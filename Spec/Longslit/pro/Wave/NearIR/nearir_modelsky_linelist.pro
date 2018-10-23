;+
; NAME:
;   nearir_modelsky_linelist
;
; PURPOSE:
;   Genreate a model sky in the near-IR with OH and H20 lines
;
; CALLING SEQUENCE:
;
; INPUTS:
;   resolution -  Instrument resolution (including slit) 
;
; OPTIONAL INPUTS:
;  T_BB -- Temperature of the BB [default=250K]
;  SCL_BB -- Scale the Blackbody [default=1]
;  SCL_OH -- Scale the OH lines
;  /NOWRITE -- Writes the FITS file but not the linelist
;  /NO_CONTI_BB -- Do not use the BB spectrum as the continuum [NOT RECOMMENDED]
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   ;
; REVISION HISTORY:
;   
;-
;------------------------------------------------------------------------------
;; nearir_modelsky_linelist, 2700., 'TSPEC_linelist.lst',
;; 'TSPEC_modelsky.fits', dlam=41.1, flgd=1
;; nearir_modelsky_linelist, 10000., 'XSHOOTER_linelist.lst',
;; 'XSHOOTER_modelsky.fits', dlam=10.0, flgd=1, wvmnx = [0.8,2.6]
;; nearir_modelsky_linelist, 540, 'GNIRS_linelist.lst',
;; 'GNIRS_modelsky.fits', dlam=100.0, flgd=1, wvmnx = [0.8,2.6]

PRO NEARIR_MODELSKY_LINELIST, resolution, linefile, outfile, WVMNX=wvmnx $
                              , INSTRUMENT = instrument, DLAM = dlam, FLGD = flgd $
                              , NOWRITE = nowrite, SCL_BB = scl_bb, T_BB = t_bb $
                              , BBONLY = bbonly $
                              , SCL_OH = scl_oh, SCL_H2O = scl_h2o, NO_CONTI_BB = no_conti_bb $
                              , pkwdth = pkwdth1, toler = toler1, thin = thin1 $
                              , fweight = fweight1, nsig = nsig1, ICLSE = ICLSE1 $
                              , WAVE_WATER = WAVE_WATER

  
  if not keyword_set(WVMNX) then wvmnx = [1., 2.6] ;; microns
  if not keyword_set(INSTRUMENT) then instrument = ''
  if not keyword_set(SCL_BB) then scl_bb = 1.
  if not keyword_set(SCL_OH) then scl_oh = 1.
  if not keyword_set(SCL_H2O) then scl_h2o = 10.
  if not keyword_set(T_BB) then T_BB  = 250.
  IF NOT KEYWORD_SET(WAVE_WATER) THEN WAVE_WATER = 2.3d

  wv_min = wvmnx[0]
  wv_max = wvmnx[1]

  ;; Peak finding stuff
  IF n_elements(pkwdth1) GT 0 THEN pkwdth = pkwdth1 ELSE pkwdth = 2.5
  IF n_elements(toler1) GT 0 THEN toler = toler1 ELSE toler = 2.5
  IF n_elements(THIN1) GT 0 THEN THIN = THIN1 ELSE THIN = 1
  IF n_elements(FWEIGHT1) GT 0 THEN FWEIGHT = FWEIGHT1 ELSE FWEIGHT = 0
  IF n_elements(ICLSE1) GT 0 THEN ICLSE = ICLSE1 ELSE ICLSE = 4L
  IF n_elements(NSIG1) GT 0 THEN NSIG = NSIG1 ELSE NSIG = 5.0d


  ;; Create the wavelength array
  case flgd of 
     0: begin ;; Linear
        ngrid = (wv_max - wv_min)/dlam
        wave = wv_min + dindgen(ngrid)*dlam 
     end
     1: begin ;; Log
        velpix = dlam  ;; km/s
        loglam = alog10(1.0d + velpix/299792.0d)
        ngrid = alog10(wv_max / wv_min) / loglam
        wave = 10^(alog10(wv_min)+findgen(ngrid)*loglam)
     end
     else: stop
  endcase
  trans = transparency(wave)
   
   ;; Empirical match to gemini broadband continuum level
   logy = - 0.55 - 0.55 * (wave-1.0)
   y = 10^logy
   
   ;; Blackbody
   planck = 6.63e-27
   c = 2.99e10
   k = 1.38e-16
   lam = wave / 1e4
   radian_per_arcsec = 1. / 3600. * 3.14159/180
   
   ;; Add in a blackbody for the atmosphere
   blackbody = 2 * planck * c^2 / lam^5 / (exp(planck*c/(lam*k*T_BB))-1)
   
   bb_counts = blackbody / (planck * c / lam) * 1e-4 * $
               (radian_per_arcsec)^2 ; * (1-trans)
   bb_counts = bb_counts * SCL_BB

;   oplot, wave, alog10(y+bb_counts), color=getcolor('red')

   ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;; add in OH lines
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
           sig = fwhm / 2.35
           gd = where(abs(wave-oh_wv[i])/sig LT 3.0, ngd)
           if (ngd GT 0) then begin
              ohspec[gd] += oh_fx[i]/20000.0 * $
                             1.0/(sqrt(2.0*!dpi)*sig)*exp(-(wave[gd]-oh_wv[i])^2/(2*sig^2))
               ;exp(-(wave[gd]-oh_wv[i])^2/(2*sig^2))*(resolution/1000.0)
           endif
       endif
    endfor
   
   ;; ;;;;;;;;;;;;;;;;;;;;;;;
   ;; H20 lines
   if wvmnx[1] GT wave_water and not keyword_set(BBONLY) then begin
      h2o_file = getenv('LONGSLIT_DIR') + '/pro/Wave/NearIR/Water/HITRAN.txt'
      readcol, h2o_file, h2o_nu, h2o_rad, format='D,D'
      h2o_wv = 1./ h2o_nu * 1e4 ;; microns
      ;; Reverse
      h2o_wv = reverse(h2o_wv)
      h2o_rad = reverse(h2o_rad) * 5e11

      ;; Cut to WVMNX
      cutwv = where( h2o_wv GT (wvmnx[0]-0.1) and $
                     h2o_wv LT (wvmnx[1]+0.1), n_h2o)
      h2o_wv = h2o_wv[cutwv]
      h2o_rad = h2o_rad[cutwv]
      
      ;; 
      mn_wv = mean(wvmnx>wave_water)
      mn = min(abs(mn_wv-h2o_wv), imn)
      h2o_dwv = abs(h2o_wv[imn] - h2o_wv[imn+1])

      dwv = mn_wv/resolution/2.3548 ;; sigma
      nsmooth = dwv / h2o_dwv

      ;; Convolve to the instrument resolution.  This is only
      ;; approximate.
      kernel = gauss_kernel(nsmooth)
      smooth_fx = convol(h2o_rad, kernel)

      ;; Rebin
      x_specrebin, h2o_wv*1e4, smooth_fx, wave*1e4, h2o_spec, /flambda
      ;; Zero out below 2.3microns (reconsider)
      zero = where(wave LT wave_water, nzero)
      if nzero GT 0 then h2o_spec[zero] = 0.
   endif else h2o_spec = fltarr(ngrid) 

   clr = getcolor(/load)
   sky_model = y+bb_counts+ohspec*SCL_OH+h2o_spec*scl_h2o
   ;stop
   x_specplot, sky_model, wav = wave, inflg = 4, /block
   
   case flgd of 
      0: begin ;; Linear
         sxaddpar, hdr, "CRVAL1", wv_min
         sxaddpar, hdr, "CDELT1", dlam
         sxaddpar, hdr, "DC-FLAG", 0
      end
      1: begin ;; Log
         sxaddpar, hdr, "CRVAL1", alog10(wv_min)
         sxaddpar, hdr, "CDELT1", loglam
         sxaddpar, hdr, "DC-FLAG", 1
      end
      else: stop
   endcase
   mwrfits, sky_model, outfile, hdr, /create
   mwrfits, wave, outfile

   
   if keyword_set(NOWRITE) then return

   ;; Now fit the continuum of the arc for accurate peak finding
   case instrument of
      'ISAAC': gdi = isaac_irwave_mask(wave) ;; hard-wired clean region mask
      else: gdi = lindgen(n_elements(wave))
   end
   ngd = n_elements(gdi)

   if keyword_set(NO_CONTI_BB) then begin
      afunc = 'LEGEND'
      npoly = 8
      tmpfit = x_fitrej(wave[gdi], sky_model[gdi], afunc, NPOLY, $
                        ivar = replicate(1., ngd), NITER = 1, $
                        hsigma = 3., lsigma = 3., rms = rms, FITSTR = fitstr, $
                        /SILENT)
      autofit = x_calcfit(wave, FITSTR = fitstr)   
   endif else autofit = y + bb_counts

   ;; Find peaks in synthetic spectrum
   x_fndpeaks, sky_model, peak, NSIG = NSIG, /silent, PKWDTH = pkwdth $
               , THIN = THIN, FWEIGHT = FWEIGHT, TOLER = TOLER $
               , ICLSE = ICLSE $
               , AUTOFIT = autofit, RMS = rms
   
   pkwave = interpol(wave, dindgen(n_elements(wave)), peak)
   m_pkwave = interpol(sky_model, dindgen(n_elements(wave)), peak)
   
   x_specplot, sky_model, wav = wave, ytwo = m_pkwave, two_wave = pkwave $
               , psym2 = 1, /block

   ;; Linelist
   label = replicate(' OH', n_elements(peak))
   qual = replicate(1, n_elements(peak))   
   
   ;; Water (assumes everything past 2.3 microns is Water and nothing
   ;; before)
   pk_h2o = where(pkwave GT wave_water, npw)
   if npw GT 0 then label[pk_h2o] = ' H2O'

   ;; OLD CODE
   ;;; Add H20 lines using FIRE R=6000 line list
   ;IF BAND EQ 'K' THEN BEGIN
   ;   ;; Open line lists
   ;   line_path = GETENV('XIDL_DIR') + '/Spec/Arcs/Lists/'
   ;   linelist = line_path + 'FIRE_OH_R6000.lst'
   ;   x_arclist, linelist, lines
   ;   ih20 = WHERE(lines.wave GE 23620.0d, nh20)
   ;   IF nh20 GT 0 THEN BEGIN
   ;      pkwave = [pkwave, lines[ih20].wave/1.0d4]
   ;      qual = [qual, replicate(1, nh20)]
   ;      label = [label, ' ' + lines[ih20].NAME]
   ;      m_pkwave = [m_pkwave, replicate(5.0, nh20)]
   ;   ENDIF
   ;ENDIF 

   writecol, linefile, pkwave, qual, label, m_pkwave
   
   RETURN
end
