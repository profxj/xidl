;+ 
; NAME:
; x_fndfitval   
;   Version 1.1
;
; PURPOSE:
;    Solves for x in y = f(x) given a fitted function f(x) and y
;
; CALLING SEQUENCE:
;   
;   xsolv = x_fndfitval(val, fitstr, xval,[fit], NITER=, IPX=,
;      TOLER=, NEG=, FITSTR=)
;
; INPUTS:
;   val        - Value to match
;   fitstr     - FIT structure
;   xval       - Values where fit was pre-evaluated
;   [fit]      - Values of the fit at xval (calculated if not input)
;
; RETURNS:
;   xsolv      - x position where the fit = val
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   NORD       - Required for LEGEND
;   TOLER      - Tolerance for match [default: 10^-4]
;   IPIX       - Starting pixel in xval
;   NEG        - Proceed in the negative direction
;   NITER      - Max number of iterations [default: 50]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xsolv = x_fndfitval(0.2, fitstr, xval)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  POLY_FIT
;  POLY
;  SVDFIT
;  SVLEG
;  BSPLIN
;  GAUSS
;
	; REVISION HISTORY:
;   23-Nov-2001 Written by JXP
;   13-Feb-2002 Requires fit structure
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_fndfitval, val, fitstr, xval, fit, $
                      IPX=ipx, TOLER=toler, NEG=neg, NITER=niter, $
                      SILENT=silent

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'xval = x_fndfitval(val, fitstr, xval, [fit], func=, ffit=, NORD=,'
    print, '         IPX=, TOLER=, /NEG )  [v1.1]'
    return, -1
  endif 

;  Optional Keywords

  if not keyword_set( IPX ) then ipx = 0
  if not keyword_set( NITER ) then niter = 50
  if not keyword_set( TOLER ) then toler = 1e-4 else begin
      if TOLER LT 1.E-8 then message, 'TOLER too small! Try rescaling.'
  endelse

  ; FIT
  if not keyword_set( FIT ) then $
    fit = x_calcfit(xval, FITSTR=fitstr)
      
; Error checking

  if (IPX+1 GT n_elements(xval)) then begin
      if not keyword_set( SILENT ) then $
        print, 'x_fndfitval: IPX > n_elem(xval): Returning n_elem'
      return, double(IPX)
  endif
  if (IPX EQ 0 AND keyword_set(NEG)) then begin
      if not keyword_set( SILENT ) then $
        print, 'x_fndfitval: IPX=0 and /NEG!: Returning 0.d'
      return, 0.d
  endif

;;;;
  if not keyword_set( NEG ) then begin
      step = 1
      qfin = n_elements(xval) - 1
  endif else begin
      step = -1
      qfin = 0
  endelse

; First bracket the root

  for q = ipx, qfin, step do begin
      if(q EQ qfin) then begin
          if not keyword_set( SILENT ) then $
            print, 'x_fndfitval: Warning -- val not bracketed'
          return, -1
      endif
      if ((fit[q]-val)*(fit[q+step]-val)) LE 0 then begin
          x1 = xval[q] 
          val1 = fit[q]
          x2 = xval[q+step]
          val2 = fit[q+step]
          break
      endif
  endfor

; Now Refine to TOLER

  diff = 10.^5
  iter = 1
  while(abs(diff) GT TOLER AND ITER LT NITER) do begin
      xcen = double(x1+x2)/2.
      cval = x_calcfit( xcen, FITSTR=fitstr)
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
      
      
      
