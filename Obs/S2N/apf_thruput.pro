;+ 
; NAME:
; apf_thruput
;    Version 1.0
;
; PURPOSE:
;    Estimates the throughput of APF using empirical measurements
;
; CALLING SEQUENCE:
;  thru = apf_thruput( wave)
;
; INPUTS:
;  wave=  -- Wavelength to calculate throughputs at
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  OUTDIR=  -- Name of output directory
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  showfits
;  querydss
;
; REVISION HISTORY:
;   27-Mar-2013 Written by BPH based on hires S2N code
;   13-Jun-2013 Rewritten by BPH, new thruput measurements
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;
function apf_thruput, wave

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'thru = apf_thruput(wave) [v1.0]'
      return,-1
  endif 

  ;; define location of output directory...
  dir = getenv('XIDL_DIR')
  if dir eq '' then dir = 'idl/xidl/'

  ;; Read sensitivity file
  if not keyword_set(SENS_FIL) then $
     sens_fil = dir + 'Obs/S2N/THRU_PUT_DATA/sens_APF_jan2014.fits'
;     sens_fil = dir + 'Obs/S2N/THRU_PUT_DATA/sens_APF_may2013.fits'
;     sens_fil = dir + 'Obs/S2N/THRU_PUT_DATA/sens_APF_29aug2012.fits'

  sens = xmrdfits(sens_fil, 1, /SILENT)
  thru = interpol(sens.eff, sens.wav, wave)
; thru = thru * 1.6 ; HACK!!!!!

;  thru = -156.702 + 0.0792474*wave - 0.0000118077 *wave*wave + 5.702421e-10 * wave * wave * wave
  thru = thru / 100.

  return, thru
end

