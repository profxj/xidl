;+ 
; NAME:
; x_constants
;   Version 1.0
;
; PURPOSE:
;    Returns a named constant
;
; CALLING SEQUENCE:
;   
;   cnst = x_constants(name)
;
; INPUTS:
;   name - Constant name
;
; RETURNS:
;
; OUTPUTS:
;   cnst  - Value
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   spl = x_constants('spl') 
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   19-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function x_constants, name

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'cnst = x_constants(name) [V1.0]'
    return, -1
  endif 

  ; Here we go
  case name of 
      'spl': return, 2.99D5
      else: return, -1
  endcase
end
