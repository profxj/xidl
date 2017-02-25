;+ 
; NAME: 
; igm_calc_lox
;    Version 1.1
;
; PURPOSE:
;    Calculate l(X) given an f(N,X) over an N_HI interval
;
; CALLING SEQUENCE:
;   
; INPUTS:
;
; RETURNS:
;  l(X)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  CUMUL=  -- Return the cumulative l(X) as a function of N_HI
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   July-2011 Written by JXP
;-
;------------------------------------------------------------------------------
function igm_calc_lox, fn_strct, z, NHI_min, iNHI_max, NEVAL=neval, CUMUL=cumul, $
                       LGNHI=lgNHI

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'lX = igm_calc_lox(fn_strct, z, NHI_min, [NHI_max]) [v1.0]'
    return, -1
  endif 

  if not keyword_set(NEVAL) then neval = 10000L
  if not keyword_set(iNHI_MAX) then begin 
     NHI_MAX = 23.
     infinity=1
  endif else NHI_MAX = iNHI_MAX

  nz = n_elements(z)

  ;; Brute force (should be good to ~0.5%)
  lgNHI = NHI_min + (NHI_MAX-NHI_MIN)*dindgen(neval)/(neval-1)
  dlgN = lgNHI[1]-lgNHI[0]
  
  ;; Evaluate f(N,X)
  lgfNX = fltarr(neval,nz)
  lX = fltarr(nz)
  for ii=0L,nz-1 do $
     lgfNX[*,ii] = eval_powerfn(fn_strct, lgNHI, z[ii])

  ;; Sum
  for ii=0L,nz-1 do $
     lX[ii] = total(10.d^(lgfNX[*,ii]+lgNHI)) * dlgN * alog(10.)
  if arg_present(CUMUL) then begin
     if nz GT 1 then stop ;; Have not modified this yet
     cumul = total(10.d^(lgfNX+lgNHI), /cumul) * dlgN * alog(10.)
  endif

  ;; Infinity?
  if keyword_set(INFINITY) then begin
     NEVAL2 = 1000L
     lgNHI2 = NHI_max + (99.-NHI_MAX)*dindgen(neval2)/(neval2-1)
     dlgN = lgNHI2[1]-lgNHI2[0]
     lgfNX = fltarr(neval2,nz)
     lX2 = fltarr(nz)
     for ii=0L,nz-1 do begin
        lgfNX[*,ii] = eval_powerfn(fn_strct, lgNHI2, z[ii])
        lX2[ii] = total(10.d^(lgfNX[*,ii]+lgNHI2)) * dlgN * alog(10.)
     endfor
     ;; 
     lX = lX + lX2
  endif

  if nz EQ 1 then lX = lX[0]
  return, lX
end
