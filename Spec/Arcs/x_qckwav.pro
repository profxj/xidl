; x_qckwav
;     Version 1.1
;
; PURPOSE:
;    Calculates the wavelength values for a set of x,y pairs in a
;    given order using the slopes calculated in x_fit2darc
;
; CALLING SEQUENCE:
;   
;  ywave = x_qckwav(xstart, ystart, ordr_str[q].arc_m,
;  arc_slope=arc_slope, NITER=, SLIT_DIST=)
;
; INPUTS:
;   xstart   -  x values
;   ystart   -  y values
;   arc_m    -  Slope of arc lines evaluated at each row in order
;
; RETURNS:
;
; OUTPUTS:
;   Fills in the profile0 and profile1 tags in the Order structure
;
; OPTIONAL KEYWORDS:
;  NITER     - Number of iterations (default: 4)
;
; OPTIONAL OUTPUTS:
;  ARC_SLOPE - The slope of the arc line at each x,y pair
;  SLIT_DIST - Distance of the x,y point along the slit from the
;              center 
;
; COMMENTS:
;
; EXAMPLES:
;  ywave = x_qckwav(xstart, ystart, ordr_str[q].arc_m, arc_slope=arc_slope)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   ??-??-2003 Written by SB
;  Feb-2005 Ported to XIDL by JXP
;-
;------------------------------------------------------------------------------

function x_qckwav, xstart, ystart, arc_m, arc_slope=arc_slope, $
                  niter=niter, slit_dist=slit_dist

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'wave = x_qckwav(xstart, ystart, arc_m, ARC_SLOPE=, NITER=' + $
        'SLIT_DIST= [v1.1]'
      return,-1
  endif

  ;; Optional keywords
  if NOT keyword_set(niter) then niter=4L

  ;; Spline
  nrow = n_elements(arc_m)
  ycol = dindgen(nrow)
  slope_spline = spl_init(ycol, double(arc_m))
  
  ywave = ystart
  ;; Iterate
  for islope=1,niter do begin 
      arc_slope = spl_interp(ycol, arc_m, slope_spline, ywave) 
      ywave = ystart - arc_slope*xstart 
      oldslope = arc_slope 
  endfor

  slit_dist =  xstart * sqrt(1.0 + arc_slope^2)
  return, ywave
end

