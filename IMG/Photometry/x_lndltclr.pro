;+ 
; NAME:
; x_lndltclr   
;    Version 1.1
;
; PURPOSE:
;    Returns the Landolt color of a star
;
; CALLING SEQUENCE:
;   
;   clr = x_lndltclr(color, lndltstr) 
;
; INPUTS:
;   color - String name of the color (e.g. 'BV')
;   lndlstr - Structure of the landolt star
;
; RETURNS:
;   clr - Color of the Landolt star
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
;   clr = x_lndltclr('BR', landolt)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------
function x_lndltclr, color, landolt
;
  if  N_params() LT 2  then begin 
      print, 'Syntax - ' +$
        'clr = x_lndltclr(color, landolt) '
      return, -1
  endif 

  ;; Parse on color
  case color of 
      'UB' : return, landolt.UB
      'UV' : return, landolt.UB + landolt.BV 
      'UR' : return, landolt.UB + landolt.BV + landolt.VR
      'BV' : return, landolt.BV
      'BR' : return, landolt.BV + landolt.VR
      'BI' : return, landolt.BV + landolt.VI
      'VR' : return, landolt.VR
      'VI' : return, landolt.VI
      'RI' : return, landolt.RI
      else : begin
          print, 'x_lndltclr: ', color, ' is not allowed !'
          return, -1
      end
  endcase
  
end

