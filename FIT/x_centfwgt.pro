;+ 
; NAME:
; x_centfwgt   
;   Version 1.1
;
; PURPOSE:
;  Given a 'peak', this routine will find the center of that peak in a
;  using flux weighting.
;
; CALLING SEQUENCE:
;   
;   xcen = x_centfwgt(xval, yval, [radius])
;
; INPUTS:
;   xval       - x values of the peak
;   yval       - y values of the peak
;   [RADIUS]   - radius of data around the peak for centroiding
;
; RETURNS:
;   xcen      - Center
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /SILENT -- Turn off warnings
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
; trace_fweight
;
; REVISION HISTORY:
;   22-Jun-2005 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;
function x_centfwgt, xval, yval, radius, SILENT=silent

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'xcen = x_centfwgt(xval, yval, [radius], /SILENT  [v1.0]'
    return, -1
  endif 

;  Optional Keywords

  if not keyword_set( RADIUS ) then radius = 3.

  xcen = n_elements(yval)/2.
  for k=0,19 do $
    xcen = trace_fweight(yval, xcen, 0L, radius=radius, xerr=xerr, $
                         invvar=ivar)
  ;; Deal with endpoints
  dx = xval[1] - xval[0]
  return, xval[0] + xcen*dx
end
