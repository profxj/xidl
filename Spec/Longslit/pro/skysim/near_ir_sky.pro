;pro near_ir_sky

   datafile = 'gemini_irsky.dat'
   resolution = 6000.
;   resolution = 400.

   openr, 10, datafile

   tmp1=0.
   tmp2=0.

   n=0
   while (not EOF(10)) do begin
       readf, 10, tmp1, tmp2
       n+=1
   endwhile
   close, 10

   wave = fltarr(n)
   flux = fltarr(n)

   openr, 10, datafile
   i=0
   while (not EOF(10)) do begin
       readf, 10, tmp1, tmp2
       wave[i] = tmp1
       flux[i] = tmp2
       i+=1
   endwhile
   close, 10

   wave = wave / 1000
   flux = flux * 1000 / 1e4

   flux = flux[where(wave GT 0.8 AND wave LT 2.5)]
   wave = wave[where(wave GT 0.8 AND wave LT 2.5)]

   plot, wave, (flux), /xstyle, xrange=[0.8,2.5]

;;;;;;;;;;;;;; Now go about creating our own synthetic version ;;;;;;;;

   wv_min = 0.8
   wv_max = 2.6
   velpix = 7.5
   loglam = alog10(1.0d + velpix/299792.0d)
   ngrid = alog10(wv_max / wv_min) / loglam
   wave = 10^(alog10(wv_min)+findgen(ngrid)*loglam)

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
   ohfile = 'rousselot2000.dat'
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

   oplot, wave, (y+bb_counts+ohspec), color=getcolor('red')

   sky = y+bb_counts+ohspec

   sxaddpar, hdr, "CRVAL1", alog10(wv_min)
   sxaddpar, hdr, "CDELT1", loglam

   mwrfits, sky, "FIRE_modelsky_R6000.fits", hdr, /create
   stop
end
