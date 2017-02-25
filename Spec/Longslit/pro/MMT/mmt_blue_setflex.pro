;+
; NAME:
;   mmt_blue_setflex
;
; PURPOSE:
; This script takes a sky spectrum for the MMT-B 300 line grating 
; setup and shifts it into the Paranal sky frame. It cretes an archived 
; reference sky spectrum. The MMT spectrum should have been reduced with 
; long_reduce with the flags nohelio and noflexure set to 1. 
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
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
;   10-Mar-2005  Written by S. Burles (MIT), David Schlegel (LBL), and 
;                Joe Hennawi (UC Berkeley)
;-
;------------------------------------------------------------------------------
; read in object number ten of the object struct which is on slit 4. 
scifile = '/Volumes/scr0/PAIRS/REDUX/mmt_redux/Dec_2006/slit2.0/Science/sci-henn121406.0043.fits.gz'

struct = mrdfits(scifile, 4)
obj = struct[2]
wave_obj = obj.wave_box
; MMT wavelengths are in vacuum
loglam = alog10(wave_obj)
; Use boxcar sky because it should not have any hot pixels
sky_obj  = obj.sky_box
disp_log = abs(loglam-shift(loglam, 1))
disp_log[0] = disp_log[1]
ny = n_elements(wave_obj)

; Ignore pixels blueward of 3200 since they are noise in the data and the 
; sky spectruum does not go that blue
dloglam = djs_median(disp_log)

; bluer slit #2 has object # 5 on it
; wave_in1 = struct[4].wave_opt
; sky_obj1 = struct[4].sky_opt 
IF TAG_EXIST(obj, 'FWHM_MED') THEN FWHM = obj.FWHM_MED $
ELSE fwhm = djs_median(obj.ARC_FWHM)
disp_LRIS = (FWHM/2.35482D)
sky_ref = skyspec_paranal(loglam, disp = disp_LRIS)

diff_log = max(loglam)-min(loglam)
xvector = (loglam-min(loglam))/diff_log

; Scale the counts sky spectrum to be the same as the Paranal spectrum
; This is effectively flux calibration.
solve_poly_ratio, xvector, sky_obj, sky_ref, replicate(1.0, ny) $
                  , npoly = 1, nback = 0 $
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

;pixnorm = WHERE(wave_obj GT 3911.0D AND wave_obj LT 3918.0D)
;sky_obj_corr = sky_obj_corr/total(sky_obj_corr[pixnorm])
;sky_ref_corr = sky_ref_corr/total(sky_ref_corr[pixnorm])

; Chop to isolate the line at 3900
;sky_obj_corr = sky_obj_corr*double(sky_obj_corr GT 0.05 $
;                                   AND wave_obj GT 3600.0 AND wave_obj LT 4000.0)
;sky_ref_corr = sky_ref_corr*double(sky_ref_corr GT 0.05 $
;                                   AND wave_obj GT 3600.0 AND wave_obj LT 4000.0)

MAXFLEX = 10D
step = lindgen(2*MAXFLEX) - MAXFLEX 
ipix = where(wave_obj LT 9200.0)

corr = c_correlate(sky_obj_corr[ipix], sky_ref_corr[ipix], step, /double)
xpeak = long_find_nminima(-corr, step, nfind = 1, minsep = 1, ypeak = ypeak $
                          , npeak = npeak, errcode = errcode $
                          , width = MAXFLEX/2.0 $
                          , /doplot, xplotfit = xfit, yplotfit = yfit)
xpeak = double(xpeak)      

print, xpeak, FORMAT = '(%"Flexure Shift: %5.2f pixels")'  

; now correct the wavelenghts using the air wavelengths
wave_ref_flex = obj.wave_box
wave_calib = interpol(wave_ref_flex, dindgen(ny), dindgen(ny) + xpeak)
sky_calib = sky_obj

; Archive flexure corrected sky
save, wave_calib, sky_calib $
      , file = '/Users/jhennawi/IDL/xidl/Spec/Longslit/calib/sky/mmt_sky_blue_300.sav'


; Plot things to check that it is okay
wave_vac = 10.0D^loglam
wave_vac_shift = interpol(wave_vac, dindgen(ny), dindgen(ny) + xpeak)
splot, wave_vac, sky_ref, psym = 10
soplot, wave_vac_shift, sky_obj_scale, col = 2, psym = 10



END
