;+
; NAME:
;    focus_test
;
; PURPOSE:
;    to measure Keck focus from DEIMOS TV guider images
;
; CALLING SEQUENCE:
;    focus_test,file_R,file_R1,[f_error,f_sigma, forcefwhm=,
;          minflux= ]
; 
; INPUTS:
;    file_R  -- fits file for R filter exposure, DEIMOS TV camera
;    file_R1 -- fits file for R1 filter
;    
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;    forcefwhm --  use an assigned fwhm for star detection
;    minflux   --  use only stars above a given flux
;
; OUTPUTS:
;    f_error  -- focus error (millimeters)
;    f_sigma  -- uncertainly in this measurement
;
; EXAMPLES:
;
; COMMENTS:
;    this simple test compares position of objects in the 2 filters,
;    seeking shifts in DEIMOS Y as a function of Y.  The R filter is 
;    non-vignetting, but the R1 filter obscures a portion of the
;    primary, the amount depending on Y position.  If in focus, the
;    shift in Y between R and R1 frames will be constant in Y.  If out
;    of focus, the delta Y will be Y dependent, with a nearly constant
;    tilt of DY vs Y.  The result is calibrated from data taken 05jul02.
;    Best results will require 1x1 binning of the TV, unless there are
;    many stars in the field.
;
; REVISION HISTORY:
;   md 14jul02
;   md 9sep02 -- attempts to make more robust
;----------------------------------------------------------------------
pro focus_test,  framefiler, framefiler1, f_error, f_sigma, $
          forcefwhm=forcefwhm, minflux=minflux

; tbd -make use of file numbers alone, rather than full names

  if n_elements(minflux) eq 0 then minflux = -1E5
  if n_elements(forcefwhm) eq 0 then forcefwhm = -1


  imager = readfits(framefiler, /silent)
  imager1 = readfits(framefiler1, /silent)
  
  nrows = (size(imager, /dimens))[0]
  fwhm = float(ceil(5.*nrows/1024.))  
  gsize = (size(imager, /dimen))[0]
  g_scale = .2076*(1024./gsize)   ;arcsec/pixel on guider
  
  if forcefwhm gt 0 then begin
;       guidefind, image,  xsout, ysout, xout, yout, flux=flux, $
;               fwhm=forcefwhm, outfwhm=fout ;find sources
       guidefind, imager,  ggr, fwhm=forcefwhm, outfwhm=foutr ;find sources
       guidefind, imager1, ggr1, fwhm=forcefwhm, outfwhm=foutr1 ;find sources
       
  endif else begin
       guidefind, imager,  ggr, fwhm=fwhm, outfwhm=foutr ;find sources
;       guidefind, imager,  ggr, fwhm=1.5*median(foutr), outfwhm=foutr
       guidefind, imager1, ggr1, fwhm=fwhm, outfwhm=foutr1 ;find sources
;       guidefind, imager1, ggr1, fwhm=1.5*median(foutr1), outfwhm=foutr1
  endelse

  good = intarr(n_elements(ggr))
  for i=0, n_elements(ggr)-1 do begin
    delta = g_scale^2*(ggr.xout -ggr[i].xout)^2  + (ggr.yout -ggr[i].yout)^2
    close = where(delta lt 2.5^2, nclose) ;object will have itself as neighbor
    if nclose lt 2 then good[i] = 1 ; no close pairs for this object  
  endfor
  ggr = ggr[where(good eq 1)] ; keep only isolated objects


     goodflux = where(ggr1.flux ge minflux, goodct)
     if goodct gt 0 then ggr1 = ggr1[goodflux]  $
     else begin
        message, "NO OBJECTS FULFILL CRITERIA in R1!!", /INFORM
        return
     ENDELSE

     goodflux = where(ggr.flux ge minflux, goodct)
     if goodct gt 0 then ggr = ggr[goodflux]  $
     else begin
        message, "NO OBJECTS FULFILL CRITERIA in R!!", /INFORM
        return
     ENDELSE

  goodcen = where(ggr.xsout ne ggr.xout OR ggr.ysout ne ggr.yout, goodct)
;seems backward to me?? MD
  IF goodct eq 0 then begin
        message, "NO OBJECTS FULFILL CRITERIA in R!!", /INFORM
        return
  ENDIF
  ggr = ggr[goodcen]

  goodcen = where(ggr1.xsout ne ggr1.xout OR ggr1.ysout ne ggr1.yout, goodct)
;seems backward to me?? MD
  IF goodct eq 0 then begin
        message, "NO OBJECTS FULFILL CRITERIA in R1!!", /INFORM
        return
  ENDIF
  ggr1 = ggr1[goodcen]


  print, 'Median flux: ', median(ggr.flux)
  print, 'Median fwhm: ', median(foutr)

  pairs = 0
  for i=0, n_elements(ggr1) -1 do begin
    delta = g_scale^2*(ggr.xout -ggr1[i].xout)^2  + (ggr.yout -ggr1[i].yout)^2
    close = where(delta lt 2.^2, nclose)
    if nclose gt 0 then $  ;add to list of close pairs
      pairs = (n_elements(pairs) eq 1) ? [ggr1[i].yout, ggr[close[0]].yout, $
        ggr1[i].xout, ggr[close[0]].xout ] : $
        [[pairs], [ggr1[i].yout, ggr[close[0]].yout, ggr1[i].xout, $  
          ggr[close[0]].xout ] ];4 vector for pairs
  endfor

;tbd-- needs robust fit to slope and error
;   slope = 0.

  xerror = .1 + .2*reform(pairs[0, *])/nrows ;rough error, gets worse where R1 is faint
  fit = svdfit(reform(pairs[0, *]), reform(pairs[1, *]-pairs[0, *]), 2, $
            /legendre, sigma=sigma, yfit=yfit, $
            measure_errors = xerror, chisq= chisqr  )

;  print, 'slope and error: (x1000) ', 1000.*fit[1], 1000.*sigma[1]
  print, 'number pts, chisqr of fit/d.f.: ', n_elements(pairs[0, *]), $
          chisqr/n_elements(pairs[0, *])

  f_error = -fit[1]/7.e-3   ;calibration based on 2x2 binned data analysis
  f_sigma = sigma[1]/7.e-3

  print, 'focus error (mm), sigma: ', f_error, f_sigma
  print, 'NOTE: positive focus error means focus should be reduced'

  window, 0
  plot, pairs[0, *], pairs[1, *]-pairs[0, *], xtitle='y position on TV', $
          ytitle='DY', psym=2
  oplot, pairs[0, *], yfit ;plot best fitting line
  
;test code here!!
  window, 1
  plot, [0, nrows], [0, nrows], /nodata, xtit='X', ytit='Y'
  dx = pairs[2, *] - pairs[3, *]
  dy = pairs[1, *] - pairs[0, *] ;delta x,y
;  print, 'mean dx,dy: ', mean(dx), mean(dy)
  dx = dx -mean(dx) ;0.368  ;mean(dx) for in-focus image
  dy = dy -mean(dy) ;0.113  ;mean(dy) 
;remove mean values to take out effect of filter tilts, mean at best
;focus (rough estimate based on ~45 stars on data of 5jul02
  arrow, pairs[3, *], pairs[1, *], pairs[3, *]+100.*dx, pairs[1, *]+ $
        100.*dy, /data

return
end


















