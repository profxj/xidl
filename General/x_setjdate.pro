;+ 
; NAME:
; x_setjdate   
;  Version 1.1
;
; PURPOSE:
;    Converts a string date into Julian date 
;
; CALLING SEQUENCE:
;   
;   jdate = x_setjdate(date, [UTTIME])
;
; INPUTS:
;   date     - String name of the date
;      Allowed formats:  DDmmmYY  (e.g. 22Oct99)
;                        DD-MM-YY (e.g. 20-12-99)
;                        YYYY-MM-DD (e.g. 1999-12-20)
;                        DDMMMYYYY (e.g. 20JAN2004)
;   [UTTIME] -- To create JD with hr, min, sec (e.g. 12:11:29)
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
;   julday  (IDL package)
;
; REVISION HISTORY:
;   31-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------
function x_setjdate, date, uttime
; Check
  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'jdate = x_setjdate(date, [uttime]) v(1.1)'
    return, -1
  endif 

; Optional Keyword
  if keyword_set( UTTIME ) then begin
      ipos = strpos(uttime,':')
      if ipos LT 0 then stop
      prs = strsplit(uttime,':',/extract)
      hr = prs[0]
      min = prs[1]
      sec = prs[2]
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
      9 : begin
          ; Fix day
          day = strmid(date, 0, 2) 
          mon = strmid(date, 2, 3)
          mon = x_getmonth(mon, /REVERSE)
          year = strmid(date, 5, 4)
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
