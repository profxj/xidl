;+ 
; NAME:
; esi_slitnm   
;     Version 1.0
;
; PURPOSE:
;    Process an image (bias subtract + flatten)
;      WARNING! Assumes 1 bias and 1 flat for all images
;
; CALLING SEQUENCE:
;   
;  nm = esi_slitnm( slit )
;
; INPUTS:
;   slit   -  Slit size
;
; RETURNS:
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
;   nm = esi_slitnm( slit ) 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function esi_slitnm, slit

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'nm = esi_slitnm(slit) [v1.0]'
      return, -1
  endif 
  
  case slit of 
      0.00: c_s = '00'
      0.30: c_s = '30'
      0.50: c_s = '50'
      0.75: c_s = '75'
      1.00: c_s = '10'
      1.25: c_s = '12'
      1.50: c_s = '15'
      6.00: c_s = '60'
      9.99: c_s = 'MH'
      else: stop
  endcase

  return, c_s
end

