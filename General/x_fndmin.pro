;+ 
; NAME:
; x_fndmin
;   Version 1.0
;
; PURPOSE:
;    Finds data point of an array nearest the input value
;      and returns the array member
;
; CALLING SEQUENCE:
;   
;   indx = x_fndmin(val, xdat, XVAL=)
;
; INPUTS:
;   xdat - Data
;   val  - value
;
; RETURNS:
;   indx  - Index in the array
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   XVAL - Value of the array at that index
;
; COMMENTS:
;
; EXAMPLES:
;   indx = x_fndmin( 1.0, findgen(1000))
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   19-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function x_fndmin, val, xdat, XVAL=xval

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'indx = x_fndmin(val, xdat, XVAL=)  [V1.0]'
    return, -1
  endif 

; Do it

  diff = abs(val - xdat)
  xval = min(diff, indx)
  delvarx, diff
  return, indx

end

