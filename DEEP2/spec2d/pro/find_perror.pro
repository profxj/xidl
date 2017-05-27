;+
; NAME:
;   find_perror
;
; PURPOSE:
;    top level routine to process one guide frame, find pointing error
;
; CALLING SEQUENCE:
;    find_perror,framefile,maskfile,offset,[forcefwhm=, minflux=, mirror=]
;
; INPUTS:
;    framefile -- guider image FITS file (1024 pixels preferred)
;                 (either full path or call from within current directory)
;    masknumber  -- string number of mask definition (db file) 
;                  containing mask definition information
;                  presumed to come from db directory, e.g '2101'
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;   forcefwhm=forcefwhm -- set to a positive value to force the use of
;                          a FWHM of forcefwhm pixels in guidefind
;   minflux=minflux -- set to throw out all objects with less than
;                      minflux counts
;   mirror  -- set this keyword to analyse only object in mirror region
;
; OUTPUTS:
;   offset -- DEIMOS (x,y) displacement (arcsec) to add to guider to match 
;       catalog positions                 
;
;
; RESTRICTIONS:
;
; COMMENTS:
;   This routine is intended for alignment purposes, to give error in
;   RA, dec, PA of telescope pointing.  This version searches all 3 unknowns
;   It is intended for use when rough alignment to of order an arcesec 
;   has already been achieved.
;
; REVISION HISTORY:
;   MD 12jun02, 
;    02jul 2002 -generalize files for use in analysis
;----------------------------------------------------------------------
pro find_perror,  framefile, masknumber, offset, mirror=mirror, forcefwhm=forcefwhm, minflux=minflux

   maskdir = '~/deep/targetselection/masks/'
;  framefile = '/home/marc/DEIMOS/2002jun04/cam0004.fits'
;  maskfile = '/home/marc/work/1HSmask.db.1145t.fits'
  maskfile = maskdir + '1HSmask.db.' + masknumber + '.fits'

  if n_elements(minflux) eq 0 then minflux = -1E5
  if n_elements(forcefwhm) eq 0 then forcefwhm = -1


  Xoffset = -7.25 ;arcsec  offsets between TV guider and DEIMOS X,Y 
  Yoffset = -2.5  ;arcsec   (possibly needs further refinement)

  image = readfits(framefile, headimage, /silent)
  nrows = (size(image, /dimens))[0]
  fwhm = float(ceil(5.*nrows/1024.))  
  gsize = (size(image, /dimen))[0]
  g_scale = .2076*(1024./gsize)   ;arcsec/pixel on guider
  
  if forcefwhm gt 0 then begin
;       guidefind, image,  xsout, ysout, xout, yout, flux=flux, $
;               fwhm=forcefwhm, outfwhm=fout ;find sources
       guidefind, image,  gg, fwhm=forcefwhm, outfwhm=fout ;find sources
       
  endif else begin
       guidefind, image,  gg, fwhm=fwhm, outfwhm=fout ;find sources
       guidefind, image,  gg, fwhm=1.5*median(fout), outfwhm=fout, /useatv
  endelse


     if keyword_set(mirror) then   gg = $  ;keep only upper region
          gg(where(gg.yout gt 580*(nrows/1024.))) ;needs to be more elegant

     goodflux = where(gg.flux ge minflux, goodct)
     if goodct gt 0 then gg = gg[goodflux]  $
     else begin
        message, "NO OBJECTS FULFILL CRITERIA!!", /INFORM
        return
     ENDELSE

  goodcen = where(gg.xsout ne gg.xout OR gg.ysout ne gg.yout, goodct)
;seems backward to me?? MD
;
;   lim = 0.5 ;limit for match of the two methods
;   goodcen = where(abs(gg.xsout -gg.xout) lt lim and abs(gg.ysout - gg.yout) $
;       lt lim, goodct)

  IF goodct eq 0 then begin
        message, "NO OBJECTS FULFILL CRITERIA!!", /INFORM
        return
  ENDIF
  gg = gg[goodcen]

  print, 'Median flux: ', median(gg.flux)
  print, 'Median fwhm: ', median(fout)

;  wset, 1
;  plot, xout, yout, psym=4, xr=[0, 1024.], yr=[0, 1024]
;  keep = where(flux gt 100.) ;cut on flux, arbitrary level for the moment
;  xout = xout[keep]
;  yout = yout[keep]

  gg.xout = float(gsize)-gg.xout ;flip X coords 
  tilt = atan((302.7-298.4)/209.)
  cos_tilt = cos(tilt)
  sin_tilt = sin(tilt)

;take out the modest rotation of guider field (rotate about guider=[0,0])
  rxout =  (gg.xout)*cos_tilt + (gg.yout)*sin_tilt 
  ryout = -(gg.xout)*sin_tilt + (gg.yout)*cos_tilt 
  DXtv = rxout*g_scale +1. ;DEIMOS X
  DYtv = -ryout*g_scale +298.4 ; DEIMOS Y (note sign inversion)

  DXtv = DXtv + Xoffset
  DYtv = DYtv + Yoffset
  
;Hi Marc,
;        Here are the corners of the TV field, in arcsec in the focal plane
;(our mask centers are 0,270 in this system):
;
;Mirror:
;  -1   94
;  -1  174
; 208  174
; 208   94

;Mask:
;  -1  174
;  -1  298.4
; 208  302.7
; 208  174
;
;Despite the decimal values on some numbers, these are probably good to about
;1 pixel
;guider plate scale: .2076 arcsec/pixel


  maskinfo = mrdfits(maskfile, 1, /silent) ;mask info
  objects  = mrdfits(maskfile, 2, /silent) ;image catalog
;retrieve bright objects from guide box region
  guideobjs = objects[where(objects.objclass eq 'G')] ;and objects.mag lt 21.)]
  
  dg_mag = guideobjs.mag
  dg_ra  = guideobjs.ra_obj - maskinfo.ra_fld ;delta ra of guide obj.
  dg_dec = guideobjs.dec_obj -maskinfo.dec_fld ;delta dec of guide. obj.

  dg_ra = dg_ra*3600.*cos(maskinfo.dec_fld * !dtor) ;arcsec. on sky
  dg_dec = dg_dec*3600.

  cost = cos(maskinfo.pa_pnt * !dtor)
  sint = sin(maskinfo.pa_pnt * !dtor) ;intended PA of mask

  DX = dg_dec*cost + dg_ra*sint ;DEIMOS X position
  DY = -(-dg_dec*sint + dg_ra*cost) +270. ;DEIMOS Y (opposite from sky)
;N.B. using convention that 4.5' from optical center is pointing center

  wset, 0
  plot, dxtv,dytv,psym=4,symsize=1.5, xr=[-10,220],yr=[310,90]
  oplot,dx,dy,psym=2
  oplot, [-1, 208], [174, 174]

;  offset1 = discrete_correlate2d([[dxtv], [dytv]],[[dx], [dy]], $
;          lagrange=[-20., 20.], nbest=fix(n_elements(dxtv)*.5), $
;            step=2, tmax=.2, tstep=.1)
;  print, '1st pass: best offset in DEIMOS (x,y,t): ', offset1
; 
;          lagrange = [-5+offset1[0], 5+offset1[1]]
;          lagrange = lagrange(sort(lagrange))
;  offset2 = discrete_correlate2d([[dxtv], [dytv]],[[dx], [dy]], $
;          lagrange=lagrange, $
;           nbest=fix(n_elements(dxtv)*.5), step=0.5, tmax=.1, tstep=.05)
;  print, '2nd pass: best offset in DEIMOS (x,y,t): ', offset2
;
;  offset = offset2
;
;  matches = [0., 0., 0., 0.]


  if goodct gt 0 then begin
     dxtv = dxtv[goodcen]
     dytv = dytv[goodcen]
     flux = gg.flux[goodcen]
  endif

  matches = [0, 0, 0, 0] 
  for i=0, n_elements(dxtv) -1 do begin
    delta = (dx-dxtv[i])^2  + (dy-dytv[i])^2
    close = where(delta lt 2.^2, nclose)
;    mindelt = min(delta, j)
;     if (mindelt lt 2.^2) then matches = [[matches], [dxtv[i], dytv[i], $
;           dx[j], dy[j]]] ;generate matching pairs for close matches
    if nclose gt 0 then begin
      maxmag = min(dg_mag[close], ind) ;select most luminous of close matches
      matches = [[matches], [dxtv[i], dytv[i], dx[close[ind]], dy[close[ind]]]]
    endif
  endfor
  if n_elements(matches) gt 4 then begin
    npairs = (size(matches, /dimens))[1]
    matches = matches[*, 1:npairs-1] ;strip off 1st placeholder
    matches = matches- gsize/2  ;get into center of TV guider frame
    npairs = npairs-1 ;cut the number by the first
    if npairs gt 1 then begin
      sums = total(matches, 2)
      deltax = (sums[2]-sums[0])/npairs
      deltay = (sums[3]-sums[1])/npairs
      print, 'Number of pairs: ', npairs
      print
      print, 'Mean offset of closely aligned stars (arcsec): '
      print, '', deltax, deltay, form='(A20,2f10.3)'
      print, 'Median offsets: ', median(matches( 2, *)-matches(0, *)), $
         median(matches(3, *)-matches(1, *)), form='(A20,2f10.3)'
      print, 'std. deviations/std. errors: ', stdev(matches( 2, *)- $
         matches(0, *)), stdev(matches(3, *)-matches(1, *)), $ 
         stdev(matches( 2, *)-matches(0, *))/sqrt(npairs), $
         stdev(matches(3, *)-matches(1, *))/sqrt(npairs), form='(A30,4f8.3)'
      thetax = total((matches[2]-matches[0])*matches[1])/total(matches[1]^2)
      thetay =-total((matches[3]-matches[1])*matches[0])/total(matches[0]^2)
      print
      print, 'Rotational shift (degrees): ', thetax*!radeg, thetay*!radeg

  endif else print, 'insufficient matches for detailed fit'
    endif

return
end










