;+ 
; NAME:
; x_fndspln
;   Version 1.0
;
; PURPOSE:
;    Solves for x in y = f(x) given a splined function f(x) and y
;
; CALLING SEQUENCE:
;   
;   xsolv = x_fndsplin(sx,sy,val,[splin] IPX=, TOLER=, /NEG=)
;
; INPUTS:
;   sx         - Values where spline was pre-evaluated
;   sy         - Values of the spline at sx
;   val        - Value to match
;   [splin]    - Spline (calculated if necessary)
;
; RETURNS:
;   xsolv      - x position where the fit = val
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   TOLER      - Tolerance for match (default: 10^-4)
;   IPIX       - Starting pixel in sx
;   /NEG       - Proceed in the negative direction from IPIX
;   /SILENT    - Turn off warning messages
;   NITER      - Max number of iterations [default: 50L]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xsolv = x_fndspln(sx, sy, val)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  SPLIN
;
; REVISION HISTORY:
;   02-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_fndspln, sx, sy, val, splin, IPX=ipx, TOLER=toler, NEG=neg, $
                      SILENT=silent, NITER=niter

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'xval = x_fndsplin(sx, sy, val, splin, '
    print, '         IPX=, TOLER=, /NEG, /SILENT)  [v1.0]'
    return, -1
  endif 

;  Optional Keywords

  if not keyword_set( IPX ) then ipx = 0
  if not keyword_set( NITER ) then niter = 50L
  if not keyword_set( TOLER ) then toler = 1d-4 else begin
      if TOLER LT 1.d-8 then message, 'TOLER too small! Try rescaling.'
  endelse

  ; SPLINT
  if not keyword_set( SPLIN ) then splin = spl_init(sx, sy, /double)
      
; Error checking

  if (IPX+1 GT n_elements(sx)) then begin
      if not keyword_set( SILENT ) then $
        print, 'x_fndsplin: IPX > n_elem(sx): Returning n_elem'
      return, double(n_elements(sx))
  endif
;;;;

; First bracket the root

  prod = (sy-val)*(shift(sy,-1)-val)
  a = where(prod LE 0, na)
  if na EQ 0 then begin
      if not keyword_set( SILENT ) then $
        print, 'x_fndsplin: Warning -- val not bracketed'
      return, -1
  endif
  if not keyword_set( NEG ) then begin
      gda = where(a GE ipx, na)
  endif else begin
      gda = where(a LE ipx, na)
  endelse

  if na EQ 0 then begin
      if not keyword_set( SILENT ) then $
        print, 'x_fndsplin: Warning -- val not bracketed'
      return, -1
  endif
  q = a[gda[0]]
  if q EQ n_elements(sx)-1 then begin
      x1 = sx[q-1] 
      val1 = sy[q-1]
      x2 = sx[q]
      val2 = sy[q]
  endif else begin
      x1 = sx[q] 
      val1 = sy[q]
      x2 = sx[q+1]
      val2 = sy[q+1]
  endelse

; Now Refine to TOLER

  diff = 1e5
  iter = 1L
  while(abs(diff) GT TOLER AND ITER LT NITER) do begin
      xcen = double(x1+x2)/2.
      cval = spl_interp(sx, sy, splin, xcen, /double)
      diff = (cval - val)/val

      if((val-val1)*(val-cval)) LE 0 then begin
          x2 = xcen
          val2 = cval
      endif else begin
          x1 = xcen
          val1 = cval
      endelse
      iter = iter+1
  endwhile

  return, xcen

end
      
      
      
