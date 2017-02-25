;+ 
; NAME:
; x_calctwilight   
;    Version 1.0
;
; PURPOSE:
;  Calculate twilight in LST given a lat+long+date
;
; CALLING SEQUENCE:
;  
; INPUTS:
;   date     - String name of the date of the night of observing at
;              the observatory (not UT)
;      Allowed formats:  DDmmmYY  (e.g. 22Oct99)
;                        DD-MM-YY (e.g. 20-12-99)
;                        YYYY-MM-DD (e.g. 1999-12-20)
;                        DDMMMYYYY (e.g. 20JAN2004)
;   obs_str   - Observatory structure (long, latitute, time zone)
;
; RETURNS: twlight in LST (decimal hours)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   stmid -- LST at midnight (end of the input date)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_calclst
;
; REVISION HISTORY:
;  Based on code in skycalendar v5 written by John Thorstensen
;  Dartmouth College
;   10-Dec-2009 Written by JXP
;-
;------------------------------------------------------------------------------

function x_calctwilight, date, obs_str, TWILIGHT_ALT=twilight_alt, STMID=stmid

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'twilight = x_calctwilight(date, obs_str, TWILIGHT_ALT) [v1.0]'
      return, -1
  endif 

  if not keyword_set(TWILIGHT_ALT) then TWILIGHT_ALT =  -18.  ;; 18 deg twilight
  ;; Grab longitude and latitude 
  longit = obs_str.longitude/15.d 
  lat = double(obs_str.latitude)

  ;; Get mjd
  jd = x_setjdate(date)  ;; Corresponds to Noon

  ;; LST midnight (ignoring daylight savings)
  jdmid = double(jd) + obs_str.tz/24.d - 0.5d ;  /* corresponding to ut */
  stmid = x_getlst(jdmid, longit)  ;; Decimal hours
  
  ;; Get RA/DEC of the Sun (decimal hours)
  x_sunradec, jdmid, rasun, decsun

  hatwilight = x_haalt(decsun,lat,TWILIGHT_ALT) 

  ;; Evening twilight
  offset = rasun+hatwilight-stmid
  if abs(offset) GT 12 then begin
      if offset GT 0 then sgn = -1 else sgn = 1
      ntf = long( (abs(offset)-12)/24 )
      offset = offset + (ntf+1)*sgn*24
  endif

  jdetwilight = jdmid + offset/24. 
  jdetwilight = x_jdsunalt(TWILIGHT_ALT,jdetwilight,lat,longit) ;
  sid_eve = x_getlst(jdetwilight,longit) ;
;  print, 'sid_eve', long(sid_eve), (sid_eve-long(sid_eve))*60

  ;; Morning twilight
  offset = rasun-hatwilight-stmid
  if abs(offset) GT 12 then begin
      if offset GT 0 then sgn = -1 else sgn = 1
      ntf = long( (abs(offset)-12)/24 )
      offset = offset + (ntf+1)*sgn*24
  endif

  jdmtwilight = jdmid + offset/24. 
  jdmtwilight = x_jdsunalt(TWILIGHT_ALT,jdmtwilight,lat,longit) ;
  sid_mor = x_getlst(jdmtwilight,longit) ;
;  print, 'sid_mor', long(sid_mor), (sid_mor-long(sid_mor))*60

  return, [sid_eve, sid_mor]
end


