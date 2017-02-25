;+ 
; NAME:
; x_oplotcirc
;
; PURPOSE:
;   Overplot a circle on a current plot given center and radius
;
; CALLING SEQUENCE:
;   x_oplotcirc, rad, x0=, y0=, NPT=, _EXTRA=
;
; INPUTS:
;  rad -- Radius of the circle (plot units)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  x0 -- x-center of the circle [default: 0.]
;  y0 -- y-center of the circle [default: 0.]
;  NPT -- Number of points to discretize the circle [default: 1000L]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  x_oplotcirc, 2., x0=3., y0=4.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   22-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_oplotcirc, radius, x0=x0, y0=y0, NPT=npt, _EXTRA=keyforplot

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
      'x_oplotcirc, radius, x0=x0, y0=y0, _EXTRA=keyforplot [v1.1]'
    return
  endif 

  if not keyword_set( NPT ) then npt = 1000L
  if not keyword_set( x0 ) then x0 = 0.
  if not keyword_set( y0 ) then y0 = 0.

  x = radius*cos(2*!pi * findgen(npt) / float(npt-1)) + x0
  y = radius*sin(2*!pi * findgen(npt) / float(npt-1)) + y0
  ;print, 'x_oplotcirc', y0

  oplot, x, y, _EXTRA=keyforplot

end

