;+ 
; NAME:
; x_getmonth  
;
; PURPOSE:
;    Converts a number into a month
;
; CALLING SEQUENCE:
;   
;   month = x_setmonth(month)
;
; INPUTS:
;   month    -  Long value of the month 
;
; RETURNS:
;   month    -  String
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
;   month = x_getmonth(10)
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   04-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------
function x_getmonth, imonth

; Check

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'month = x_getmonth(imonth) v(1.0)'
    return, -1
  endif 

;

  case imonth of 
      1 : return, 'JAN'
      2 : return, 'FEB'
      3 : return, 'MAR'
      4 : return, 'APR'
      5 : return, 'MAY'
      6 : return, 'JUN'
      7 : return, 'JUL'
      8 : return, 'AUG'
      9 : return, 'SEP'
     10 : return, 'OCT'
     11 : return, 'NOV'
     12 : return, 'DEC'
      else : begin
          print, 'x_getmonth: Not prepared for this string', imonth
          return, -1
      end
  endcase

end
