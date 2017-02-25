;+ 
; NAME:
; x_getlst   
;    Version 1.0
;
; PURPOSE:
;  Calculate twilight in LST given a lat+long+date
;
; CALLING SEQUENCE:
;  
; INPUTS:
;   jdin     - Julian Date
;      Allowed formats:  DDmmmYY  (e.g. 22Oct99)
;                        DD-MM-YY (e.g. 20-12-99)
;                        YYYY-MM-DD (e.g. 1999-12-20)
;                        DDMMMYYYY (e.g. 20JAN2004)
;   longit   - Longitude on Earth [decimal hours]
;
; RETURNS: twlight in LST  strarr(2)
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
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_getlst
;
; REVISION HISTORY:
;  Based on code in skycalendar v5 written by John Thorstensen
;  Dartmouth College
;   10-Dec-2009 Written by JXP
;-
;------------------------------------------------------------------------------

function x_getlst, jdin, longit

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'twilight = x_getlst(mjd or DATE) [v1.0]'
      return, -1
  endif 

  ;; Get mjd
  if size(jdin,/type) EQ 7 then jd = x_setjdate(jdin) else jd = double(jdin)

  ;; LST
	;; returns the local MEAN sidereal time (dec hrs) at julian date jd
	;; at west longitude long (decimal hours).  Follows
        ;; definitions in 1992 Astronomical Almanac, pp. B7 and L2. 
        ;; Expression for GMST at 0h ut referenced to Aoki et al, A&A 105,
	;; p.359, 1982. 

  jdnoon2000jan1 = 2451545.d ;

  jdint = long(jd)
  jdfrac = jd - jdint
  if jdfrac LT 0.5 then begin
      jdmid = jdint - 0.5d       ;
      ut = jdfrac + 0.5d         ;
  endif else begin
      jdmid = jdint + 0.5d       ;
      ut = jdfrac - 0.5d         ;
  endelse

  t = (jdmid - jdnoon2000jan1)/36525 ;
  sid_g = (24110.54841+8640184.812866*t+0.093104*t*t-6.2e-6*t*t*t)/86400. ;
  sid_int = long(sid_g)               ;
  sid_g = sid_g - double(sid_int) ;
  sid_g = sid_g + 1.0027379093 * ut - longit/24. ;
  sid_int = long(sid_g)               ;
  sid_g = (sid_g - double(sid_int)) * 24. ;
  if sid_g LT 0. then sid_g = sid_g + 24. ;

  return, sid_g

end


