;+ 
; NAME:
; x_setddate   
;
; PURPOSE:
;    Converts a string date into Julian date 
;
; CALLING SEQUENCE:
;   
;   jdate = x_setjdate(date)
;
; INPUTS:
;   date     - String name of the date
;
; RETURNS:
;   jdate    -  Decimal date
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   If Year is < 50, then the year 2000 is presumed
;
; EXAMPLES:
;   jdate = x_setjdate('29Oct00')
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   31-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------
function x_setjdate, date, uttime


; Check

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'jdate = x_setjdate(date, [uttime]) v(1.0)'
    return, -1
  endif 

; Optional Keyword
  if keyword_set( UTTIME ) then begin
      hr = strmid(uttime,0,2)
      min = strmid(uttime,3,2)
      sec = strmid(uttime,6)
  endif

  slen = strlen(date)

  case slen of 
      7 : begin
          day = strmid(date, 0, 2)
          mon = strmid(date, 2, 3)
          year = strmid(date, 5, 2)
          iyr = fix(year)
          if iyr LT 50 then year = strjoin(['20',year])
      end
      8 : begin
          day = strmid(date, 0, 2)
          mon = strmid(date, 3, 2)
          year = strmid(date, 6, 2)
          iyr = fix(year)
          if iyr LT 50 then year = strjoin(['20',year])
      end
      10 : begin
          ; Fix day
          if strmid(date,8,1) EQ ' ' then day = '0'+strmid(date, 9, 1) $
          else day = strmid(date, 8, 2)
          mon = strmid(date, 5, 2)
          year = strmid(date, 0, 4)
      end
      else : begin
          print, 'x_setjdate: Not prepared for this string', date
          return, -1
      end
  endcase
  
  ; Use JULDAY
  if keyword_set( UTTIME ) then $
    jdate = julday(mon, day, year, hr, min, sec) $
  else jdate = julday(mon, day, year)

  return, jdate

end
