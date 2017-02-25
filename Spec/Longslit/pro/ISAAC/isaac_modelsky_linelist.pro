PRO ISAAC_MODELSKY_LINELIST, band1, linefile = linefile, outfile = outfile
  
  band = strcompress(band1, /rem)
  arcpath = GETENV('XIDL_DIR') + '/Spec/Arcs/Lists/'

  instrument = 'ISAAC'
  IF NOT KEYWORD_SET(OUTFILE) THEN outfile = $
     arcpath + instrument + '_modelsky_MR_' + band + '.fits' 
  IF NOT KEYWORD_SET(LINEFILE) THEN linefile = $
     arcpath + instrument + '_modelsky_OH_linelist_MR_' + band + '.lst'
  CASE instrument OF 
     'ISAAC': BEGIN
        CASE BAND OF
           ;; These resolution values for the 0.6" slit. Numbers
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
        slit_width = 0.6d ;; hard wiring everything for 0.6" slit for now
        plate_scale = 0.147d
        fnslit = slit_width/plate_scale
        pkwdth = 1.3d*fnslit
        toler = fnslit/2.0d > 2.0d
        THIN = 1
        FWEIGHT = 0
        NSIG = 5.0d
     END
     ELSE: message, 'other instruments not yet supported'
  ENDCASE
  
   ngrid = (wv_max - wv_min)/dlam
   wave = wv_min + dindgen(ngrid)*dlam 
   trans = transparency(wave)
   
   ;; Empirical match to gemini broadband continuum level
   logy = - 0.55 - 0.55 * (wave-1.0)
   y = 10^logy
   
   planck = 6.63e-27
   c = 2.99e10
   k = 1.38e-16
   lam = wave / 1e4
   radian_per_arcsec = 1. / 3600. * 3.14159/180
   
   ;; Add in a blackbody for the atmosphere

   blackbody = 2 * planck * c^2 / lam^5 / (exp(planck*c/(lam*k*250.))-1)
   
   bb_counts = blackbody / (planck * c / lam) * 1e-4 * $
               (radian_per_arcsec)^2 ; * (1-trans)
   
;   oplot, wave, alog10(y+bb_counts), color=getcolor('red')

   ;; add in OH lines
   
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
      ;; produces better wavelength solutions with 1.0 threshold
       if (oh_fx[i] GT 1.0d) then begin

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
   sky_model = y+bb_counts+ohspec
   ;;x_specplot, sky, wav = wave, inflg = 4, /block
   
   sxaddpar, hdr, "CRVAL1", wv_min
   sxaddpar, hdr, "CDELT1", dlam
   sxaddpar, hdr, "DC-FLAG", 0
   mwrfits, sky_model, outfile, hdr, /create 
   
   ;; Now fit the continuum of the arc for accurate peak finding
   IF KEYWORD_SET(NO_CONTI_BB) THEN BEGIN
      gdi = isaac_irwave_mask(wave) ;; hard-wired clean region mask
      ngd = n_elements(gdi)
      afunc = 'LEGEND'
      tmpfit = x_fitrej(wave[gdi], sky_model[gdi], afunc, 8, $
                        ivar = replicate(1., ngd), NITER = 1, $
                        hsigma = 3., lsigma = 3., rms = rms, FITSTR = fitstr, $
                        /SILENT)
      autofit = x_calcfit(wave, FITSTR = fitstr)   
   ENDIF ELSE autofit = y + bb_counts
   ;; Find peaks in synthetic OH spectrum
   x_fndpeaks, sky_model, peak, NSIG = NSIG, /silent, PKWDTH = pkwdth $
               , THIN = THIN, FWEIGHT = FWEIGHT, TOLER = TOLER $
               , AUTOFIT = autofit, RMS = rms
   
   pkwave = interpol(wave, dindgen(n_elements(wave)), peak)
   m_pkwave = interpol(sky_model, dindgen(n_elements(wave)), peak)
   x_specplot, sky_model, wav = wave, ytwo = m_pkwave, two_wave = pkwave $
               , psym2 = 1, /block
   label = replicate(' OH', n_elements(peak))
   qual = replicate(1, n_elements(peak))   
   ;; Add H20 lines using FIRE R=6000 line list
   IF BAND EQ 'K' THEN BEGIN
      ;; Open line lists
      line_path = GETENV('XIDL_DIR') + '/Spec/Arcs/Lists/'
      linelist = line_path + 'FIRE_OH_R6000.lst'
      x_arclist, linelist, lines
      ih20 = WHERE(lines.wave GE 23620.0d, nh20)
      IF nh20 GT 0 THEN BEGIN
         pkwave = [pkwave, lines[ih20].wave/1.0d4]
         qual = [qual, replicate(1, nh20)]
         label = [label, ' ' + lines[ih20].NAME]
         m_pkwave = [m_pkwave, replicate(5.0, nh20)]
      ENDIF
   ENDIF 

   writecol, linefile, pkwave, qual, label, m_pkwave
   
   RETURN
end
