;+ 
; NAME:
; x_mode   
;    Version 1.1
;
; PURPOSE:
;  Calcualtes the mode after one inputs an array (default: nearest integer)
;
; CALLING SEQUENCE:
;   
; 
;
; INPUTS:
;   arr  -- Array of values
;
; RETURNS:
;   mode
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
;   mode = x_mode(array)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   25-Aug-2004 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_mode, arr

  ; 
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mode = x_mode(arr) [v1.1]'
      return, -1
  endif 

  hist = histogram(arr)
  mx = max(hist, imx)
  return, float(imx)

end
