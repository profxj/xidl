;+ 
; NAME:
; x_addnoise
;   Version 1.0
;
; PURPOSE:
;  Simply bins data up in integer pixels    
;
; CALLING SEQUENCE:
;   
;   bin = x_addnoise(fx, nbin)
;
; INPUTS:
;   fx       - Flux
;   nbin     - Number of pixels to bin on
;
; RETURNS:
;   bin       - Structure of data
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   bin = x_addnoise(fx, 3)
;
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
