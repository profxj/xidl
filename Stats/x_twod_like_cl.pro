;+ 
; NAME:
; x_twod_like_cl
;    Version 1.0
;
; PURPOSE:
;  Find the confindence levels for an input 2D likelihood function
;  Assumes a fixed, uniform grid
;
; CALLING SEQUENCE:
;  DeltaL = x_twod_like_cl( lnL, [cl], SIGMA=sigma)
;
; INPUTS:
;   lnL :: Evaluated log likelihood (normalized to zero)
;   [cl] :: Confidence level to consider
;
; RETURNS:
;  
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  SIGMA=  Can be input in lieu of CL
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Dec-2004 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_twod_like_cl, lnL, CL, SIGMA=sigma, GAUSSV=gaussv

  ; 
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'err = x_twod_like_cl( lnL, [CL], SIGMA=sigma) [v1.0]'
      return, 0.
  endif 
  if not keyword_set(CL) then begin
     if not keyword_set(SIGMA) then return,-1
     cl = 1.- 2*(1.-gauss_pdf(sigma))
  endif
     
  ;; Normalize and sort
  mx = max(lnL)
  norm_lnL = lnL - mx
  srt = sort(lnL)
  norm_lnL = norm_lnL[srt]

  ;; Area under the surface
  cumul_area = total( exp(norm_lnL > (-15.d)), /cumul)
  tot_area = max(cumul_area)
  cumul_area = cumul_area/tot_area

  DeltaL =  -1*interpol(norm_lnL, cumul_area, 1-cl, /spline)

  ;; Calcualte the Gaussian one?
  if arg_present(gaussv) then begin
     if n_elements(cl) EQ 1 then $
        gaussv = chisqr_cvf(1.-cl, 2) / 2.
  endif

  return, DeltaL

end
