;+
; NAME:
;   long_setflex
;
; PURPOSE: This script takes a sky spectrum and shifts it into frame
; of an archived sky spectrum taken from Paranal. It creates an
; archived reference sky spectrum of AIR wavelengths versus sky, which
; is written to outfile. The input scifile spectrum should have been
; reduced with long_reduce with the flags nohelio=1 and noflexure=1.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;  /LICKSKY -- Use the Lick sky instead of Paranal
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
;
; REVISION HISTORY:
;   17-Aug-2007  Written Joe Hennawi (UC Berkeley)
;-
;------------------------------------------------------------------------------
pro long_setflex, scifile, outfil, sciind = sciind, npoly = npoly $
                  , LICKSKY = licksky $
                  , WAVE_MIN = WAVE_MIN, WAVE_MAX = WAVE_MAX


  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'long_setflex, scifil, outfil [v1.0]'
      return
  endif 

  IF NOT KEYWORD_SET(NPOLY) THEN NPOLY = 5
  IF NOT KEYWORD_SET(sciind) THEN sciind = 1L
  IF NOT KEYWORD_SET(WAVE_MIN) THEN WAVE_MIN = 3100.0d
  IF NOT KEYWORD_SET(WAVE_MAX) THEN WAVE_MAX = 10000.0d

  hdr = xheadfits(scifile)
  instrument = strtrim(sxpar(hdr, 'INSTRUME'))
  grism = strtrim(sxpar(hdr, 'GRISNAME'))
  grating = strtrim(sxpar(hdr, 'GRANAME'))
  IF instrument EQ 'LRISBLUE' AND grism EQ '1200/3400' THEN LRISB1200 = 1
  IF instrument EQ 'LRIS' AND grating EQ '1200/7500' THEN LRISR1200 = 1
  struct = xmrdfits(scifile, 5)
  obj = struct[sciind-1L]
  wave_obj = obj.wave_box
; wavelengths are in vacuum
  loglam = alog10(wave_obj)
; Use boxcar sky because it should not have any hot pixels
  sky_obj  = obj.sky_box
  disp_log = abs(loglam-shift(loglam, 1))
  disp_log[0] = disp_log[1]
  ny = n_elements(wave_obj)

  dloglam = djs_median(disp_log)

  IF TAG_EXIST(obj, 'ARC_FWHM_MED') THEN FWHM = obj.ARC_FWHM_MED $
  ELSE fwhm = djs_median(obj.ARC_FWHM_FIT)
  disp_LRIS = (FWHM/2.35482D)
  if keyword_set(LICKSKY) then $
    sky_ref = skyspec_lick(loglam, disp = disp_LRIS) $
  else sky_ref = skyspec_paranal(loglam, disp = disp_LRIS)
  
  diff_log = max(loglam)-min(loglam)
  xvector = (loglam-min(loglam))/diff_log

  
; Scale the counts sky spectrum to be the same as the Paranal spectrum
; This is effectively flux calibration.
  solve_poly_ratio, xvector, sky_obj, sky_ref, replicate(1.0, ny) $
                    , npoly = npoly, nback = 0 $
                    , yfit = sky_obj_scale, ymult = ymult2, yadd = yadd2
; Subtract off smooth spectrum from each
  obj_set = bspline_iterfit(wave_obj, sky_ref $
                            , nord = 3, upper = 3.0, lower = 3.0 $
                            , nbkpts = 10, yfit = sky_ref_cont)
  obj_set = bspline_iterfit(wave_obj, sky_obj_scale $
                            , nord = 3, upper = 3.0, lower = 3.0 $
                            , nbkpts = 10, yfit = sky_obj_scale_cont)
  sky_obj_corr = sky_obj_scale - sky_obj_scale_cont
  sky_ref_corr = sky_ref - sky_ref_cont
  
  ;; Use these next few lines for the LRIS-B with the 1200 line grism
  IF KEYWORD_SET(LRISB1200) THEN BEGIN
      pixnorm = WHERE(wave_obj GT 3911.0D AND wave_obj LT 3918.0D)
      sky_obj_corr = sky_obj_corr/total(sky_obj_corr[pixnorm])
      sky_ref_corr = sky_ref_corr/total(sky_ref_corr[pixnorm])
      
      djs_iterstat, sky_obj_corr, mean = mean, sigma = sigma
      chop = mean + 6.0*sigma
      ;; Chop to isolate the line at 3900. 
      sky_obj_corr = sky_obj_corr*double(sky_obj_corr GT chop AND $
                                         wave_obj GT 3600.0 AND  $
                                         wave_obj LT 4000.0)
      sky_ref_corr = sky_ref_corr*double(sky_ref_corr GT chop AND $
                                         wave_obj GT 3600.0 AND  $
                                         wave_obj LT 4000.0)
  ENDIF ELSE IF KEYWORD_SET(LRISR1200) THEN BEGIN
      pixnorm = WHERE(wave_obj GT 5574.0D AND wave_obj LT 5586.0D)
      sky_obj_corr = sky_obj_corr/total(sky_obj_corr[pixnorm])
      sky_ref_corr = sky_ref_corr/total(sky_ref_corr[pixnorm])
;      djs_iterstat, sky_obj_corr, mean = mean, sigma = sigma
;      chop = mean + 2.0*sigma
;      sky_obj_corr = sky_obj_corr*double(sky_obj_corr GT chop AND $
;                                         ((wave_obj GT 5574.0 AND  $
;                                           wave_obj LT 5586.0) OR  $
;                                          (wave_obj GT 5886.0 AND  $
;                                           wave_obj LT 5902.0)))
;      sky_ref_corr = sky_ref_corr*double(sky_ref_corr GT chop AND $
;                                         ((wave_obj GT 5574.0 AND  $
;                                           wave_obj LT 5586.0) OR  $
;                                          (wave_obj GT 5886.0 AND  $
;                                           wave_obj LT 5902.0)))
  ENDIF
  MAXFLEX = 60D
  step = lindgen(2*MAXFLEX) - MAXFLEX 
  ipix = where(wave_obj GT WAVE_MIN AND wave_obj LT WAVE_MAX)
  corr = c_correlate(sky_obj_corr[ipix], sky_ref_corr[ipix], step, /double)
  xpeak = long_find_nminima(-corr, step, nfind = 1, minsep = 1, ypeak = ypeak $
                            , npeak = npeak, errcode = errcode $
                            , width = MAXFLEX/2.0 $
                            , /doplot, xplotfit = xfit, yplotfit = yfit)
  print, xpeak, FORMAT = '(%"Flexure Shift: %5.2f pixels")'  
  
; now correct the wavelenghts using this shift
  wave_ref_flex = obj.wave_box
  wave_calib = interpol(wave_ref_flex, dindgen(ny), dindgen(ny) + xpeak)
  sky_calib = sky_obj

; Archive flexure corrected sky
  save, wave_calib, sky_calib, file =outfil

; Plot things to check that it is okay
  wave_vac = 10.0D^loglam
  wave_vac_shift = interpol(wave_vac, dindgen(ny), dindgen(ny) + xpeak)
  x_specplot, sky_ref, wav = wave_vac, inflg = 4 $
              , ytwo = sky_obj_scale, two_wave = wave_vac_shift $
              , title = 'Black is archive spectrum. Purple is shifted sky data' 
  ;splot, wave_vac, sky_ref, psym = 10
  ;soplot, wave_vac_shift, sky_obj_scale, col = 2, psym = 10

  return
END
