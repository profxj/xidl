;+ 
; NAME:
; x_haalt   
;    Version 1.0
;
; PURPOSE:
;	Returns hour angle at which object at dec is at altitude alt */
;
; CALLING SEQUENCE:
;  
; INPUTS:
;   dec    - deg
;   lat    - Latitude on Earth [deg]
;   alt    - Altitude relative to the horizon 
;
; RETURNS: HA of the object
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
;  x_haalt
;
; REVISION HISTORY:
;  Based on code in skycalendar v5 written by John Thorstensen
;  Dartmouth College
;   10-Dec-2009 Written by JXP
;-
;------------------------------------------------------------------------------
function x_haalt, idec, ilat, alt

  ;; Define
  DEG_IN_RADIAN = 180.d/!dpi
  HRS_IN_RADIAN = DEG_IN_RADIAN/15.d
  PI = !dpi
  dec = idec
  lat = ilat

;	min_max_alt(lat,dec,&min,&max);
;	if(alt < min) 
;		return(1000.);  /* flag value - always higher than asked */
;	if(alt > max)
;		return(-1000.); /* flag for object always lower than asked */
  dec = (0.5*PI) - dec / DEG_IN_RADIAN ;
  lat = (0.5*PI) - lat / DEG_IN_RADIAN ;
  coalt = (0.5*PI) - alt / DEG_IN_RADIAN ;
  x = (cos(coalt) - cos(dec)*cos(lat)) / (sin(dec)*sin(lat)) ;
  if abs(x) LE 1. then return, acos(x) * HRS_IN_RADIAN $;
  else stop ;printf("Error in ha_alt ... acos(>1).\n") ;
  
  return, -1
end

;#define PI                3.14159265358979
;#define ARCSEC_IN_RADIAN  206264.8062
;#define DEG_IN_RADIAN     57.2957795130823
;#define HRS_IN_RADIAN     3.819718634
;#define J2000             2451545.
;#define SEC_IN_DAY        86400.
;#define FLATTEN           0.003352813 
;#define EQUAT_RAD         6378137.
;#define TWILIGHT_ALT      -18.
