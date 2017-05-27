;+
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; splinefwhm2.pro
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; PURPOSE
;     The function splinefwhm.pro takes a radial profile (e.g. flux
;     as a function of radius) and uses a spline fit to determine the
;     full-width at half-maximum (fwhm) of the profile. The radial
;     profile is interpolated at 50x as many points so as to more
;     easily fit the spline function. Then using the fit, the fwhm is
;     taken to occur at the smallest radius at which the spline-fit
;     profile drops below half of the maximum. 
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; SYNTAX
;     fwhm = splinefwhm(rad, prof, splrad, splprof)
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; INPUTS
;     rad = a vector conatining the radii.
;     prof = a vector containing the counts at each radius given in
;            rad. 
;     splrad = a variable which will be set equal to the radii at
;              which the spline function was fit. 
;     splprof = a variable which will be set equal to the fit spline
;               profile. 
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; KEYWORDS
;     None.
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; OUTPUTS
;     fwhm = the derived full-width at half-maximum value according to
;            the spline fit. A failure in the fit will return a fwhm
;            value of 999. 
;     splrad, splprof (see INPUTS)
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; PROCEDURES CALLED 
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; EXAMPLES
;     None.
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; HISTORY
;     Created August 1, 1995 by M. Liu.
;     Revised June 24, 2002 by mcc.
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;-

FUNCTION splinefwhm2, rad, prof, splrad, splprof

;;; CHECK THAT THE NUMBER OF PARAMETERS IS > 2.
IF n_params() LT 2 THEN BEGIN
	PRINT,'function splinefwhm2, radius, profile, [splrad, splprof]'
	RETALL
ENDIF

;;; CHECK THAT rad AND prof HAVE THE SAME NUMBER 
;;; OF ELEMENTS. 
nrad = N_ELEMENTS(rad)
IF nrad NE N_ELEMENTS(prof) THEN $
   PRINT,'radius and profile have unequal # of elements'

;;; IF THE rad OR prof VECTORS HAVE < 3 ELEMENTS, 
;;; THEN RETURN fwhm=999.
IF N_ELEMENTS(rad) LT 3 THEN BEGIN
   fwhm = 999.
   RETURN, fwhm
ENDIF

;;; INTERPOLATE THE RADIAL COORDINATES (rad) AT 50
;;; TIMES AS MANY POINTS TO PROVIDE ADDITIONAL POINTS
;;; AT WHICH TO FIT THE SPLINE.
splrad = MIN(rad) + FINDGEN(nrad*50+1) * (MAX(rad)-MIN(rad)) / (nrad*50)
nspl = N_ELEMENTS(splrad)

;;; USE THE IDL PROCEDURE SPLINE TO FIT A SPLINE 
;;; PROFILE TO THE RADIAL PROFILE. 
splprof = SPLINE(rad,prof,splrad)

;;; FIND WHERE THE SPLINE PROFILE DROPS TO 1/2 
;;; OF THE MAXIMUM. IF THERE ARE NO POINTS IN THE
;;; SPLINE PROFILE LESS THAN THE HALF-MAXIMUM, 
;;; THEN SET THE INDEX i EQUAL TO THE NUMBER 
;;; OF ELEMENTS IN THE VECTOR rad. ELSE SET
;;; THE INDEX i EQUAL TO THE LOWEST VALUE 
;;; IN THE INDEX VECTOR wh. 
wh = WHERE(splprof LT 0.5*MAX(splprof), whct)
IF whct EQ 0 THEN i = nspl ELSE i = MIN(wh)

;;; IF WE DID NOT FIND A POINT LESS THAN 
;;; HALF-MAXIMUM OR WE FOUND IT TOO CLOSE 
;;; TO THE MAXIMUM RADIUS (1 AWAY), THEN
;;; RETURN THE VALUE fwhm=999.
IF (i LT 2) OR (i EQ nspl) THEN BEGIN
  RETURN, 999.
ENDIF

;;; SUM THE 2 POINTS STRADDLING THE HALF-MAX 
;;; THEREBY ESTIMATING THE WIDTH (fwhm).
fwhm = splrad(i) + splrad(i-1)

RETURN, fwhm

END



