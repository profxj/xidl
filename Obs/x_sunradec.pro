;+ 
; NAME:
; x_sunradec   
;    Version 1.0
;
; PURPOSE:
;  
;
; CALLING SEQUENCE:
;  
; INPUTS:
;
; RETURNS: 
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
;  x_calclst
;
; REVISION HISTORY:
;  Based on code in skycalendar v5 written by John Thorstensen
;  Dartmouth College
;   10-Dec-2009 Written by JXP
;-
;------------------------------------------------------------------------------

;/* Low precision formulae for the sun, from Almanac p. C24 (1990) */
;/* ra and dec are returned as decimal hours and decimal degrees. */

pro x_sunradec, jd, ra, dec

  DEG_IN_RADIAN = 180.d/!dpi
  HRS_IN_RADIAN = DEG_IN_RADIAN/15.d
  PI = !dpi

  n = jd - 2451545.0d             ;

  L = 280.460d + 0.9856474 * n   ;
  g = (357.528d + 0.9856003 * n)/DEG_IN_RADIAN ;
  lambda = (L + 1.915 * sin(g) + 0.020 * sin(2. * g))/DEG_IN_RADIAN ;
  epsilon = (23.439 - 0.0000004 * n)/DEG_IN_RADIAN ;

  x = cos(lambda)               ; 
  y = cos(epsilon) * sin(lambda) ; 
  z = sin(epsilon)*sin(lambda)  ;

  ;; Calculate theta
  if x EQ 0  then begin
      if y GT 0. then theta = PI/2.  $ ;
      else if y LT 0. then  theta = 3.*PI/2. $
           else theta = 0.      ;   /* x and y zero */
  endif else theta = atan(y/x);
  if x LT 0. then theta = theta + PI ;
  if theta LT 0. then theta = theta + 2.*PI ;

  ra = theta*HRS_IN_RADIAN ;
  dec = (asin(z))*DEG_IN_RADIAN ;

  return
end
