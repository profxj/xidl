;+ 
; NAME: 
; x_tefflya_norm   
;    Version 1.1
;
; PURPOSE:
;    Normalize an f(N,z) distribution with teff Lya
;
; CALLING SEQUENCE:
;   
;  x_initpowerfn, powerfn
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;  powerfn -- Structure containing a power-law f(N,X)
;
; OPTIONAL KEYWORDS:
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
;
; REVISION HISTORY:
;   May-2011 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Initalize the powerfn structure using empirical measures
pro x_tefflya_norm, powerfn_strct, idx, tau0, EW_spline

  if  N_params() LT 4  then begin 
     print,'Syntax - ' + $
           'x_tefflya_norm, powerfn_strct, idx, tau0, EW_spline [v1.0]'  
     return
  endif 
  if not keyword_set(NHI_MIN) then NHI_MIN = 11.5
  if not keyword_set(NHI_MAX) then NHI_MAX = 22.0 
  if not keyword_set(N_eval) then N_eval = 5000L

  ;; ;;;;
  ;; Set initial normalization
  powerfn_strct[idx].fn_pivot[0] = -99.
  powerfn_strct[idx].fn_pivot[1] = 0.  ;; Log

  for ii=1L,powerfn_strct[idx].npivot do begin
     powerfn_strct[idx].fn_pivot[ii+1] = powerfn_strct[idx].fn_pivot[ii] + $
                                         (powerfn_strct[idx].pivots[ii+1] - $
                                          powerfn_strct[idx].pivots[ii])* $
                                         powerfn_strct[idx].beta[ii] 
  endfor

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; "Integrate" the 'B' parameter

  lgNval = NHI_MIN + (NHI_MAX-NHI_MIN)*dindgen(N_eval)/(N_eval-1) ;; Base 10 
  dlgN = lgNval[1]-lgNval[0]
  Nval = 10.d^lgNval

  ;; Get log f(N,z) but without z evolution
  tmpfn = powerfn_strct
  tmpfn[*].gamma[*] = 0.
  logfNz = eval_powerfn(tmpfn, lgNval, 0., FNSTR=tmpfn[idx])

  restEW = spl_interp(EW_SPLINE.NHI, EW_SPLINE.EW, EW_SPLINE.splint, lgNval) ;; Ang

  Bprim = total( 10.d^(logfNz + lgNval + alog10(restEW))) * dlgN * alog(10.)  / 1215.6701
  scale = alog10(tau0[0]/Bprim)


  powerfn_strct[idx].fn_pivot[1:powerfn_strct[idx].npivot+1] = $
     powerfn_strct[idx].fn_pivot[1:powerfn_strct[idx].npivot+1] + scale 
  powerfn_strct[idx].fn_pivot[powerfn_strct[idx].npivot+2:*] = -99.
;  printcol, powerfn_strct[idx].pivots, powerfn_strct[idx].fn_pivot, $
;            powerfn_strct[idx].beta
;  stop
  
  return
end
