;+ 
; NAME:
; x_addnoise
;   Version 1.1
;
; PURPOSE:
;  Simple routine to add noise to a normalized spectrum.  
;
; CALLING SEQUENCE:
;   
;   newf = x_addnoise(fx, snr, SEED=)
;
; INPUTS:
;   fx       - Flux array
;   snr      - Signal-to-noise per pixel of the data
;
; RETURNS:
;   newfx    - Flux array with noise
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  SEED=  -- Seed for random number generator [Default: -1322]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   newfx = x_addnoise(fx, 15.)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Nov-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_addnoise, fx, snr, SEED=seed

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'newfx = x_addnoise(fx, snr) [V1.0]'
    return, -1
  endif 

  if not keyword_set(SEED) then seed = -1322L

; Create noise array

  ;; Double?
  if size(fx, /type) EQ 5 then flg_dbl = 1 else flg_dbl = 0
  
  npix = n_elements(fx)
  return, fx + randomn(seed, npix, DOUBLE=flg_dbl, /normal)/snr

end
