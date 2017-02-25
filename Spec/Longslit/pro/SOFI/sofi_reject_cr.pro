;PRO gnirs_reject_cr, a1, ivar_a1, a2, ivar_a2, tset_slits, piximg, ordermask $
;                     , cr_mask1 = cr_mask1, cr_mask2 = cr_mask2 $
;                     , CR_THRESH = CR_THRESH, CR_VAL = CR_VAL

path = '/Users/joe/DATA/SOFI_DATA/2011-09-22/'
TELLURIC = 0
;init = 175
init = 278
nseq = 2
IF KEYWORD_SET(TELLURIC) AND nseq NE 1 THEN message, 'Do not stack Tellurics'
even = 2*lindgen(nseq)
odd = 2*lindgen(nseq) + 1
anum = init + even
bnum = init + odd
afiles = path + 'SOFI_' + string(anum, FORMAT = '(I4.4)') + '.fits'
bfiles = path + 'SOFI_' + string(bnum, FORMAT = '(I4.4)') + '.fits'

slit_arcsec = 0.60d
plate_scale = 0.288d
slit = slit_arcsec/plate_scale        
FWHM = slit
sigma_psf = FWHM/2.35482D
pkwdth = slit
TOLER = slit/2.0D
nseq = n_elements(afiles)
;fnseq = float(nseq)
;fnseq2 = fnseq*fnseq

FOR ii = 0L, nseq-1L DO BEGIN
   sofi_proc, afiles[ii], aimg, ivar_a, hdr = hdr_a
   sofi_proc, bfiles[ii], bimg, ivar_b, hdr = hdr_b
   IF ii EQ 0 THEN BEGIN
      dims = size(aimg, /dim)
      nx = dims[0]
      ny = dims[1]
      a_stk = fltarr(nx*ny, nseq)
      ivar_a_stk = fltarr(nx*ny, nseq)
      mask_a_stk = fltarr(nx*ny, nseq)
      b_stk = fltarr(nx*ny, nseq)
      ivar_b_stk = fltarr(nx*ny, nseq)
      mask_b_stk = fltarr(nx*ny, nseq)
      hdr = hdr_a
   ENDIF
   a_stk[*, ii] = aimg
   ivar_a_stk[*, ii] = ivar_a
   mask_a_stk[*, ii] = (ivar_a LE 0.0)
   b_stk[*, ii] = bimg
   ivar_b_stk[*, ii] = ivar_b
   mask_b_stk[*, ii] = (ivar_b LE 0.0)
ENDFOR
;; Create an average sky for the wavelengths
skystack = fltarr(nx*ny, 2*nseq)
maskstack = fltarr(nx*ny, 2*nseq)
skystack[*, 0:(nseq-1L)] = a_stk[*, 0:nseq-1L]
skystack[*, nseq:(2*nseq-1L)] = b_stk[*, 0:nseq-1L]
maskstack[*, 0:(nseq-1L)] =  mask_a_stk[*, 0:nseq-1L]
maskstack[*, nseq:(2*nseq-1L)] = mask_b_stk[*, 0:nseq-1L]

avg_sky =  reform(djs_avsigclip(skystack, 2, sigrej = sigrej, inmask = maskstack, outmask = outmask) $
                  , nx, ny)
smashmask = reform((total(outmask, 2) EQ 2*nseq), nx, ny) 
finalmask = (smashmask EQ 0)

tset_slits = niri_slitset(nx, ny)
slitmask = long_slits2mask(tset_slits)
ximg = long_slits2x(tset_slits)
y_img = replicate(1.0, nx)#findgen(ny)

IF NOT KEYWORD_SET(BSP) THEN BSP = 0.8
IF NOT KEYWORD_SET(CR_THRESH) THEN CR_THRESH = 10.0
IF NOT KEYWORD_SET(CR_VAL) THEN CR_VAL    = 1000.0

pixset = long_wavepix(avg_sky, tset_slits, FWHM = FWHM, pkwdth = pkwdth, toler = toler, chk = chk_trc)
piximg = long_wpix2image(pixset, tset_slits)

fnnow = float(nseq-1L)
fnnow2 = fnnow*fnnow
ivec = lindgen(nseq)
FOR inow = 0L, nseq-1L DO BEGIN
   ;; Create an average image and associated noise from all but this image
   ind_now = ivec[WHERE(ivec NE inow)]
   aden = total((mask_a_stk[*, ind_now] EQ 0), 2) 
   a_bar = fnnow*(aden GT 0.0)*total(a_stk[*, ind_now]*(mask_a_stk[*, ind_now] EQ 0), 2)/ $
           (aden + (aden EQ 0.0))
   var_a_stk = (ivar_a_stk 0.0)/(ivar_a_stk + (ivar_a_stk LE 0.0))
   var_a_bar = fnnow2*(aden GT 0.0)*total(var_a_stk[*, ind_now]*(mask_a_stk[*, ind_now]EQ 0), 2)/ $
               (aden + (aden EQ 0.0))^2 
   ivar_a_bar = (var_a_bar GT 0.0)/(var_a_bar + (var_a_bar LE 0.0))


   adiff_cr = 



ENDFOR
stop

sky_resids = fltarr(nx, ny)
y_img = findgen(ny)## replicate(1.0, nx)

adiff_cr = (a1-a2)
var_a1   = (ivar_a1 GT 0.0)/(ivar_a1 + (ivar_a1 LE 0.0))
var_a2   = (ivar_a2 GT 0.0)/(ivar_a2 + (ivar_a2 LE 0.0))
var_adiff  = var_a1 + var_a2
ivar_adiff = (var_adiff GT 0.0)/(var_adiff + (var_adiff LE 0.0))

FOR iorder = 0L, norders-1L DO BEGIN
    inorder = where(ordermask EQ order_vec[iorder], nord)
    fitpix = where(ordermask EQ order_vec[iorder] AND $
                   finite(adiff_cr)  AND $
                   edgmask EQ 0 AND $
                   abs(adiff_cr) LE 5.0d4 AND $
                   ivar_adiff GT 0.0, npix)
    psort = sort(piximg[fitpix])
    sset = bspline_iterfit(piximg[fitpix[psort]], adiff_cr[fitpix[psort]] $
                           , invvar = (ivar_adiff[fitpix[psort]] GT 0.0) $
                           , upper = 3, lower = 3 $
                           , bkspace = bsp, maxiter = 20, maxrej = 10 $
                           , /silent)
    sky_resids[inorder] = bspline_valu(piximg[inorder], sset)
ENDFOR
adiff_final = adiff_cr - sky_resids
;   Take the PSF width to be that of the spectra direction (3 pixels) for 
;   the 3-pixel slit. This prevents the routine from rejecting sky lines
;   This is a description of the 3x3 core of the 2D PSF for reject_cr.pro
;    
;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
;    PSFVALS[0]          1.   PSFVALS[0]
;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
FWHM = 3.0                  
sigma_psf = FWHM/2.35482D
psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
cr_sharp1 = psf_reject_cr(adiff_final, ivar_a1, psfvals $
                           , satmask = (adiff_final LT -10))
cr_sharp2 = psf_reject_cr(-adiff_final, ivar_a2, psfvals $
                          , satmask = (-adiff_final LT -10))
;; chi^2 CR rejection. The CRs have to be positive, which keeps from 
;; us from masking persistence which can be negative. 
X2_adiff = adiff_cr*sqrt(ivar_adiff)
;cr_mask1 = lonarr(nx, ny)
;cr_mask2 = lonarr(nx, ny)
;FOR iorder = 0L, norders-1L DO BEGIN
;    ipix = WHERE(ordermask EQ order_vec[iorder])
;    djs_iterstat, X2_adiff[ipix], mean = mean, sigma = sigma, sigrej = 2.0
;    cr_mask1[ipix] = ((X2_adiff[ipix] GT (mean + CR_thresh*sigma)) AND $
;                      (a1[ipix] GT CR_VAL))
;    cr_mask2[ipix] = ((X2_adiff[ipix] LT (mean - CR_thresh*sigma)) AND $
;                      (a2[ipix] GT CR_VAL))
;ENDFOR
ipix = WHERE(ordermask GT 0)
djs_iterstat, X2_adiff[ipix], median = median, sigma = sigma, sigrej = 2.0
cr_mask1 = cr_sharp1 OR (X2_adiff GT  median + CR_thresh*sigma AND a1 GT CR_VAL)
cr_mask2 = cr_sharp2 OR (X2_adiff LT  median - CR_thresh*sigma AND a2 GT CR_VAL)

cr_mask1 = long(cr_mask1 EQ 0)
cr_mask2 = long(cr_mask2 EQ 0)

RETURN
END
