;+ 
; NAME:
;  grb_autochrt
;   Version 1.1
;
; PURPOSE:
;    Creates a DSS chart from a BACODINE.  [old and not well tested]
;
; CALLING SEQUENCE:
;   grb_autochrt, fil
;
; INPUTS:
;     fil -- Email containing the coord
;
; RETURNS:
;   flux= -- PS file of the finding chart
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
;   grb_fxlum, 'baco_050730'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Jul-2005 Written by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro grb_autochrt, fil, _EXTRA=extra

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'grb_autochrt, fil,  [v1.1]'
    return
  endif 

  ;; Fil
  if strlen(fil) LT 2 then begin
      a = findfile('baco*', count=count)
      case count of
          0: begin
              print, 'grb_autochrt:  Specify the file!'
              return
          end
          1: 
          else: begin
              print, 'grb_autochrt: Multiple baco found'
              print, 'grb_autochrt: Using ', a[0]
          end
      endcase
      fil = a[0]
  endif

  ;; Parse
  readcol, fil, tag, val, format='A,F'

  ;; Create finding chart
  
  ;; Date
  date = where(strtrim(tag,2) EQ 'GRB_DATE:')
  tjd = val[date]
  daycnv, tjd+2440000.5, yr, mn, day
  if mn LT 10 then mnc = '0'+strtrim(mn,2) else mnc = strtrim(mn,2)
  if day LT 10 then dayc = '0'+strtrim(day,2) else dayc = strtrim(day,2)

  nm = 'GRB'+strmid(strtrim(yr,2),2)+mnc+dayc

  ;; RA/DEC
  rat = where(strtrim(tag,2) EQ 'POINT_RA:')
  dect = where(strtrim(tag,2) EQ 'POINT_DEC:')
  x_radec, ras, decs, val[rat[0]], val[dect[0]], /flip

  x_fndchrt, [nm,ras,decs], /RADEC, _EXTRA=extra

  return
end
