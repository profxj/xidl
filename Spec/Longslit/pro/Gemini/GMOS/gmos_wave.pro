PRO gmos_wave, databasefile, idfile, gsfile, slitfile, wavefile  $
               , hdr = hdr, xshift = xshift, yshift = yshift $
               , FWHMSET = FWHMSET, ROTATION=rotation $
               , QAFILE = QAFILE

;    Transform wavelengths stored in the gemini mosaic files
;     to 3 images representing wavelengths in the raw CCD frames
;

if N_elements(xshift) NE 3 then begin
    xshift = [-2.50, 0, 3.5]
    print, 'Using these xshifts: ', xshift
endif

if N_elements(yshift) NE 3 then begin
    yshift = [-1.58, 0, -1.0]
    print, 'Using these yshifts: ', yshift
endif

if not keyword_set(ROTATION) then rotation = [-0.020, 0, 0.04] * !dPi/180.

lines = djs_readlines(databasefile)

if long(lines[6]) NE 1 then begin
    print, 'Need to implement non-chebyshev basis functions.'
    return
endif

xcoeff = long(lines[7])
ycoeff = long(lines[8])


xmin = float(lines[10]) - 1.
xmax = float(lines[11]) - 1.
ymin = float(lines[12]) - 1.
ymax = float(lines[13]) - 1.

coeff = dblarr(xcoeff, ycoeff)
for i = 14, (14+xcoeff*ycoeff-1L) DO $
  coeff[i-14] = double(lines[i])

xbin = long(2048 * 3. / (xmax + 1))     + 1
ybin = long(4600 / (ymax + 1))     + 1


ncol = 2048L/xbin
nrow = 4608L/ybin
rawwave = dblarr(ncol, nrow, 3)

;
;  let's do the middle image first, no shift or rotation
;
xstart2 = long((xmax+1) /2 - 2048/xbin/2)

x = 2*(dindgen(ncol) + xstart2 - xmin)/(xmax-xmin) - 1. 
y = 2*(dindgen(nrow) - ymin)/(ymax-ymin) - 1.

xbasis = fchebyshev(x, xcoeff)
ybasis = fchebyshev(y, ycoeff)

rawwave[*, *, 1] = (xbasis # coeff) # transpose(ybasis)

;
;  Now let's do the other CCDs
;
for i = 0, 2, 2 do begin
    
;  make starting x,y images
    
    x = (dindgen(ncol) # replicate(1, nrow))[*] * xbin - ((ncol*xbin-1)/2)
    y = (dindgen(nrow) ## replicate(1, ncol))[*] * ybin - ((nrow*ybin-1)/2)
    
    
    c = cos(rotation[i]) 
    s = sin(rotation[i]) 
    xrot = x * c + y * s
    yrot = y * c - x * s

    xbeg = 1.0d * xstart2 *xbin + xshift[i] + (i-1.) * (37 + ncol*xbin)
    xs = (xrot + ((ncol*xbin-1)/2) + xbeg)/xbin
    ys = (yrot + ((nrow*ybin-1)/2) + yshift[i])/ybin
    
    xnorm = 2*(xs - xmin)/(xmax-xmin) - 1. 
    ynorm = 2*(ys - ymin)/(ymax-ymin) - 1. 
    
    xbasis = fchebyshev(xnorm, xcoeff)
    ybasis = fchebyshev(ynorm, ycoeff)

    rawwave[*, *, i] = 0.
    for iy = 0, ycoeff-1 do $
      for ix = 0, xcoeff-1 do $
        rawwave[*,*,i] = rawwave[*,*,i] + $
      xbasis[*, ix] * ybasis[*, iy] * coeff[ix, iy]
    
endfor 
waveimg = gmos_mosaic(rawwave)
;  Now fill in the gaps with a polynomial fit
dims = size(waveimg)
nx = dims[1]
ny = dims[2]
xc = replicate(1.0, nx) # findgen(ny)
ivar = dblarr(nx, ny) +1.0
inmask = lonarr(nx, ny) + 1L
zero_inds = WHERE(waveimg LE 0.0, nzero)
ivar[zero_inds] = 0.0
inmask[zero_inds] = 0
;  Should be using chebyshev here but it breaks for some reason
xy2traceset, transpose(xc), transpose(waveimg), wset, ncoeff = 5 $
             , inmask = transpose(inmask), invvar = transpose(ivar) $
             , func = 'poly', upper = 3, lower = 3, /silent
traceset2xy, wset, transpose(xc), tfit
fit = transpose(tfit)
waveimg[zero_inds] = fit[zero_inds]

qafile = repstr(wavefile, '.fits', '.ps')
; If the gsfile and idfile are set fit for the arc shapes
stop
IF KEYWORD_SET(GSFILE) AND KEYWORD_SET(IDFILE) THEN $
  gmos_arcshapes, gsfile, idfile, waveimg, slitfile $
  , FWHMSET = FWHMSET, REGION = REGION, MAXSEP = MAXSEP, QAFILE = QAFILE

;   Construct a pixel image from the wavelength image for GEMINI
wave_max = max(waveimg[WHERE(waveimg GT 0.0)])
wave_min = min(waveimg[WHERE(waveimg GT 0.0)])
dispersion =  (wave_max-wave_min)/double(ny)
piximg = (waveimg GT 0.0)*(waveimg-wave_min)/dispersion

;mwrfits, waveimg, wavefile, hdr, /create 
mwrfits, waveimg, wavefile, /create 
mwrfits, piximg, wavefile 
mwrfits, fwhmset, wavefile 

return
end

