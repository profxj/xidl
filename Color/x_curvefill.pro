;+ 
; NAME:
; x_curvefill
;
; PURPOSE:
;    Shade the region between two curves
;
; CALLING SEQUENCE:
;      x_curvefill, x, y1, y2, COLOR=, OUTLINECOLOR=, OUTHICK=
;   
; INPUTS:
;  x -- x-values of the curve
;  y1 -- lower set of y-values
;  y2 -- upper set of y-values
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; color=  -- Color for shading
; /LINE_FILL -- Hashes instead of fills
; ORIENTATION -- Sets the angle of the lines
; /xyreverse -- Lets you fill based on x rather than y
;               ie: x->y, y1->x1, y2->x2 
; x2=  -- Lets you input a second array specifying both x1 and x2.
; TRANSPARENT= 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   18-May-2006 Written by Joe H.
;   7-Dec-2009 Modified by Marc R. (xyreverse and x2)
;-
;------------------------------------------------------------------------------
pro x_curvefill, x, y1, y2, color = color, outlinecolor = outlinecolor $
                 , _EXTRA = EXTRA, OUTHICK = OUTHICK, xyreverse=xyreverse, x2=x2

  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_curvefill, x, y1, y2, COLOR=, /LINE_FILL, ORIENTATION=, TRNASPARENT= [v1.1]'
    return
  endif 

  n1 = n_elements(x)
  IF n_elements(y1) NE n1 OR n_elements(y2) NE n1 or $
    N_elements(y1) NE n_elements(y2) THEN message, 'problem with array sizes'

  if keyword_set(xyreverse) then begin
     ;; Create a polygon to fill from x values
     xpoly = [y1[0], y2, reverse(y1)]
     ypoly = [x[0], x, reverse(x)]
  endif else begin
     ;; Create a polygon to fill from y values (original code)
     xpoly = [x[0], x, reverse(x)]
     ypoly = [y1[0], y2, reverse(y1)]
  endelse 
  
  ; This is in case you have two curves with different x1, x2, y1, y2
  ; If this is specified, it overwrites the polygons from anything before.
  if keyword_set(x2) then begin
     if n_elements(x2) ne n_elements(x) then message, 'problem with array sizes'
     x1 = x
     xpoly=[x1[0], x2, reverse(x1)]
     ypoly=[y1[0], y2, reverse(y1)]
  endif

  
  ;PolyFill, xpoly, ypoly, Color = color, NOCLIP = 0, _EXTRA = EXTRA
  cgcolorfill, xpoly, ypoly, Color = color, NOCLIP = 0, _EXTRA = EXTRA

  IF n_elements(OUTLINECOLOR) NE 0 THEN $
    plots, xpoly, ypoly, Color = outlineColor, Thick = OUTHICK, NOCLIP = 0
  return
  
END
