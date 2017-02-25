;+ 
; NAME:
; x_maxcl
;    Version 1.1
;
; PURPOSE:
;  Find the confindence levels for an input likelihood function:
;  lnL(x).  One dimension only.
;
; CALLING SEQUENCE:
;  err = x_maxcl( x, lnL, maxL, CL=cl )
;
; INPUTS:
;   x :: Range of x-values for the likelihood function
;   lnL :: Evaluated log likelihood (normalized to zero)
;   maxL :: Value of x where lnL is maximal
;
; RETURNS:
;  
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  CL=  -- Confidence limit (Default=0.683)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_maxcl
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Dec-2004 Written by JXP
;-
;------------------------------------------------------------------------------
function x_maxcl_spln, x

common x_maxclcmm, mxcl_x, mxcl_y, mxcl_sp, mxcl_cl, mxcl_x0, mxcl_flg

    return, spl_interp(mxcl_x, mxcl_y, mxcl_sp, x, /double)
end

;;;
function x_maxcl_func, x

common x_maxclcmm, mxcl_x, mxcl_y, mxcl_sp, mxcl_cl, mxcl_x0, mxcl_flg

  if mxcl_flg EQ 1 then $
    return, abs(qromb('x_maxcl_spln',mxcl_x0,x,/double)-mxcl_cl) $
  else $
    return, abs(qromb('x_maxcl_spln',x,mxcl_x0,/double)-mxcl_cl) 

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_maxcl, x, lnL, maxL, CL=cl

common x_maxclcmm, mxcl_x, mxcl_y, mxcl_sp, mxcl_cl, mxcl_x0, mxcl_flg

  ; 
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'err = x_maxcl( x, lnL, maxL, CL=) [v1.1]'
      return, 0.
  endif 
  if not keyword_set(CL) then cl = 0.683
  mxcl_cl = cl
  mxcl_x0 = maxL

  ;; Integrate first
  total = int_tabulated(x,exp(lnL),/sort,/double)

  ;; Spline
  mxcl_x = x
  mxcl_y = exp(lnL) / total
  mxcl_sp = spl_init(mxcl_x, mxcl_y, /double)

  ;; Fractions
  xmin = min(mxcl_x, max=xmax)
  f1 = qromb('x_maxcl_spln',xmin,mxcl_x0,/double)
  f2 = qromb('x_maxcl_spln',mxcl_x0,xmax,/double)
  
  err = dblarr(2)

  ;; Upper
  mxcl_flg=1
  mxcl_cl = cl*f2
  err[0] = x_golden('x_maxcl_func',mxcl_x0,(mxcl_x0+xmax)/2.,xmax)
  ;; Lower
  mxcl_flg=0
  mxcl_cl = cl*f1
  err[1] = x_golden('x_maxcl_func',xmin,(mxcl_x0+xmin)/2.,mxcl_x0)

  return, err

end
