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
;   sy         - 
;   val        - Value to match
;   [splin]      - Spline
;
; RETURNS:
;   xsolv      - x position where the fit = val
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   TOLER      - Tolerance for match (default: 10^-4)
;   IPIX       - Starting pixel in sx
;   NEG        - Proceed in the negative direction
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
  if not keyword_set( TOLER ) then toler = 10.0^(-4) else begin
      if TOLER LT 1.E-8 then message, 'TOLER too small! Try rescaling.'
  endelse

  ; SPLINT
  if not keyword_set( SPLIN ) then splin = spl_init(sx, sy, /double)
      
; Error checking

  if (IPX+1 GT n_elements(sx)) then begin
      if not keyword_set( SILENT ) then $
        print, 'x_fndsplin: IPX > n_elem(sx): Returning n_elem'
      return, double(n_elements(sx))
  endif
  if (IPX EQ 0 AND keyword_set(NEG)) then begin
      if not keyword_set( SILENT ) then $
        print, 'x_fndsplin: IPX=0 and /NEG!: Returning 0.d'
      return, 0.d
  endif

;;;;
  if not keyword_set( NEG ) then begin
      step = 1
      qfin = n_elements(sx) - 1
  endif else begin
      step = -1
      qfin = 0
  endelse

; First bracket the root

  for q = ipx, qfin, step do begin
      if(q EQ qfin) then begin
          if not keyword_set( SILENT ) then $
            print, 'x_fndsplin: Warning -- val not bracketed'
          return, -1
      endif
      if ((sy[q]-val)*(sy[q+step]-val)) LE 0 then begin
          x1 = sx[q] 
          val1 = sy[q]
          x2 = sx[q+step]
          val2 = sy[q+step]
          break
      endif
  endfor

; Now Refine to TOLER

  diff = 10.^5
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
      
      
      
