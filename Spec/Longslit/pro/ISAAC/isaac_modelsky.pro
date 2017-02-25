PRO ISAAC_MODELSKY, band, outfile = outfile
  
  IF NOT KEYWORD_SET(OUTFILE) THEN outfile = 'ISAAC_modelsky_MR_' + band + '.fits' 
  CASE BAND OF
     ;; These resolut9ion values for the 0.6" slit. Numbers
     ;; taken from ISAAC manual Table 6
     'H': BEGIN
        resolution = 5100.0d
        dlam = (0.079d/1024.0d)
        wv_min = 1.3
        wv_max = 1.9
      END
     'K': BEGIN
        resolution = 4400.0d
        dlam = (0.122d/1024.0d)
        wv_min = 1.9
         wv_max = 2.6
      END
  ENDCASE
  
 

   ngrid = (wv_max - wv_min)/dlam
   wave = wv_min + dindgen(ngrid)*dlam 
   trans = transparency(wave)
   
   ; Empirical match to gemini broadband continuum level
   logy = - 0.55 - 0.55 * (wave-1.0)
   y = 10^logy


   planck = 6.63e-27
   c = 2.99e10
   k= 1.38e-16
   lam = wave / 1e4
   radian_per_arcsec = 1. / 3600. * 3.14159/180

  ; Add in a blackbody for the atmosphere

   blackbody = 2 * planck * c^2 / lam^5 / (exp(planck*c/(lam*k*250.))-1)

   bb_counts = blackbody / (planck * c / lam) * 1e-4 * $
     (radian_per_arcsec)^2 ; * (1-trans)

;   oplot, wave, alog10(y+bb_counts), color=getcolor('red')

   ; add in OH lines

;   ohfile = '/home/simcoe/idl/GNIRS/pro/ohlines_rousselot.dat'
   ohfile = getenv('XIDL_DIR') + '/Spec/Longslit/pro/skysim/rousselot2000.dat'
   openr, 10, ohfile
   noh = 0
   while (not EOF(10)) do begin
       readf, 10, tmp1, tmp2
       noh+=1
   endwhile
   close, 10

   oh_wv = fltarr(noh)
   oh_fx = dblarr(noh)

   openr, 10, ohfile
   i=0
   while (not EOF(10)) do begin
       readf, 10, tmp1, tmp2
       oh_wv[i] = tmp1/10000.
       oh_fx[i] = 1*tmp2
       i+=1
   endwhile
   close, 10

   ohspec = fltarr(ngrid)
   for i=0, noh-1 do begin
       if (oh_fx[i] GT 1.) then begin

           fwhm = oh_wv[i]/resolution
           sigma = fwhm / 2.35
           gd = where(abs(wave-oh_wv[i]) LT 0.001, ngd)
           if (ngd GT 0) then begin
               ohspec[gd] += oh_fx[i] / 160. * $
;                 1/(sqrt(2*3.14159)*sigma) * $
                 exp(-(wave[gd]-oh_wv[i])^2/(2*sigma^2)) * (resolution/1000)
           endif
       endif
   endfor
   clr = getcolor(/load)
   sky = y+bb_counts+ohspec
   x_specplot, sky, wav = wave, inflg = 4, /block
   
   IF KEYWORD_SET(OUTFILE) THEN BEGIN
   sxaddpar, hdr, "CRVAL1", wv_min
   sxaddpar, hdr, "CDELT1", dlam
   sxaddpar, hdr, "DC-FLAG", 0


   mwrfits, sky, outfile, hdr, /create 

end
