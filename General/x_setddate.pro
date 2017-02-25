;+ 
; NAME:
; x_setddate   
;  Version 1.1
;
; PURPOSE:
;    Converts a string date into decimal date 
;
; CALLING SEQUENCE:
;   ddate = x_setddate(date)
;
; INPUTS:
;   date     - String name of the date
;      Allowed formats:  DDmmmYY  (e.g. 22Oct99)
;                        DD-MM-YY (e.g. 20-12-99)
;                        YYYY-MM-DD (e.g. 1999-12-20)
;
; RETURNS:
;   ddate    -  Decimal date
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
;   ddate = x_setddate('29Oct00')
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   07-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------
function x_setddate, date

; Check

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'ddate = x_setddate(date) v(1.1)'
    return, -1
  endif 

;

  slen = strlen(date)

  case slen of 
      7 : begin
          day = strmid(date, 0, 2)
          mon = strmid(date, 2, 3)
          year = strmid(date, 5, 2)
          iyr = fix(year)
          if iyr LT 50 then year = strjoin(['20',year])

          goddate = strjoin([day,'-',mon,'-',year])

          rgdate = date_conv(goddate, 'REAL')

          yr = fix(rgdate/1000)
          dy = (rgdate-double(yr)*1000)/356.

          ddate = double(yr) + dy
      end
      8 : begin
          day = strmid(date, 0, 2)
          mon = strmid(date, 3, 2)
	  mon = x_getmonth(long(mon))
          year = strmid(date, 6, 2)
          iyr = fix(year)
          if iyr LT 50 then year = strjoin(['20',year])

          goddate = strjoin([day,'-',mon,'-',year])

          rgdate = date_conv(goddate, 'REAL')

          yr = fix(rgdate/1000)
          dy = (rgdate-double(yr)*1000)/356.

          ddate = double(yr) + dy
      end
      9 : begin
          day = strmid(date, 0, 2)
          mon = strmid(date, 3, 2)
	  mon = x_getmonth(long(mon))
          year = strmid(date, 6, 2)
          iyr = fix(year)
          if iyr LT 50 then year = strjoin(['20',year])

          goddate = strjoin([day,'-',mon,'-',year])

          rgdate = date_conv(goddate, 'REAL')

          yr = fix(rgdate/1000)
          dy = (rgdate-double(yr)*1000)/356.

          ddate = double(yr) + dy
      end
      10 : begin
          ; Fix day
          if strmid(date,8,1) EQ ' ' then day = '0'+strmid(date, 9, 1) $
          else day = strmid(date, 8, 2)
          mon = strmid(date, 5, 2)
	  mon = x_getmonth(long(mon))
          year = strmid(date, 0, 4)
;          iyr = fix(year)
;          if iyr LT 50 then year = strjoin(['20',year])

          goddate = strjoin([day,'-',mon,'-',year])

          rgdate = date_conv(goddate, 'REAL')

          yr = fix(rgdate/1000)
          dy = (rgdate-double(yr)*1000)/356.

          ddate = double(yr) + dy
      end
      21 : begin
          ; Fix day
          if strmid(date,8,1) EQ ' ' then day = '0'+strmid(date, 9, 1) $
          else day = strmid(date, 8, 2)
          mon = strmid(date, 5, 2)
	  mon = x_getmonth(long(mon))
          year = strmid(date, 0, 4)
;          iyr = fix(year)
;          if iyr LT 50 then year = strjoin(['20',year])

          goddate = strjoin([day,'-',mon,'-',year])

          rgdate = date_conv(goddate, 'REAL')

          yr = fix(rgdate/1000)
          dy = (rgdate-double(yr)*1000)/356.

          ddate = double(yr) + dy
      end
      else : begin
          print, 'x_setddate: Not prepared for this string', date
          return, -1
      end
  endcase

  return, ddate

end
