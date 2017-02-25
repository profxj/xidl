;+ 
; NAME:
; x_getmonth  
;
; PURPOSE:
;    Converts a number into a month
;
; CALLING SEQUENCE:
;   
;   month = x_setmonth(ival)
;
; INPUTS:
;   ival    -  Long value of the month 
;
; RETURNS:
;   month    -  String
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /REVERSE  --  Turn a month into a number
;  /SSMALL  --  Use lowercase
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
function x_getmonth, imonth, SMALL=small, REVERSE=reverse, SSMALL=ssmall

; Check

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'month = x_getmonth(imonth) v(1.1)'
    return, -1
  endif 

  if keyword_set(REVERSE) then begin
      case imonth of
          'JAN': return, 1
          'Jan': return, 1
          'jan': return, 1
          'FEB': return, 2
          'Feb': return, 2
          'feb': return, 2
          'MAR': return, 3
          'Mar': return, 3
          'mar': return, 3
          'APR': return, 4
          'Apr': return, 4
          'apr': return, 4
          'MAY': return, 5
          'May': return, 5
          'may': return, 5
          'JUN': return, 6
          'Jun': return, 6
          'jun': return, 6
          'JUL': return, 7
          'Jul': return, 7
          'jul': return, 7
          'AUG': return, 8
          'Aug': return, 8
          'aug': return, 8
          'SEP': return, 9
          'Sep': return, 9
          'sep': return, 9
          'OCT': return, 10
          'Oct': return, 10
          'oct': return, 10
          'NOV': return, 11
          'Nov': return, 11
          'nov': return, 11
          'DEC': return, 12
          'Dec': return, 12
          'dec': return, 12
          else: return, -1
      endcase
  endif

;
  if not keyword_set( SMALL ) and not keyword_set(SSMALL) then begin
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
              print, 'x_getmonth: Not prepared for this value', imonth
              return, -1
          end
      endcase
  endif else begin
      if keyword_set(SMALL) then begin
          case imonth of 
              1 : return, 'Jan'
              2 : return, 'Feb'
              3 : return, 'Mar'
              4 : return, 'Apr'
              5 : return, 'May'
              6 : return, 'Jun'
              7 : return, 'Jul'
              8 : return, 'Aug'
              9 : return, 'Sep'
              10 : return, 'Oct'
              11 : return, 'Nov'
              12 : return, 'Dec'
              else : begin
                  print, 'x_getmonth: Not prepared for this value', imonth
                  return, -1
              end
          endcase
      endif 
      if keyword_set(SSMALL) then begin
          case imonth of 
              1 : return, 'jan'
              2 : return, 'feb'
              3 : return, 'mar'
              4 : return, 'apr'
              5 : return, 'may'
              6 : return, 'jun'
              7 : return, 'jul'
              8 : return, 'aug'
              9 : return, 'sep'
              10 : return, 'oct'
              11 : return, 'nov'
              12 : return, 'dec'
              else : begin
                  print, 'x_getmonth: Not prepared for this value', imonth
                  return, -1
              end
          endcase
      endif 
  endelse

end
