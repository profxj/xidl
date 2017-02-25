;+ 
; NAME:
; x_qsolumfunc
;
; PURPOSE:
;    Returns the luminosity function of quasars for a default (or
;    specified) range of luminosity or magnitude at a given redshift.
;    Several options are provided
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
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
;   Nov-2011 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_qsolumfunc, z, flg, M1450_EVAL=m1450_eval 

  if (N_params() LT 2) then begin 
     print,'Syntax - ' + $
           'phi = x_qsolumfunc(z,flg, M1450_EVAL=) [v1.0]'
     return, -1
  endif 
  
  if not keyword_set(M1450_EVAL) then $
     M1450_EVAL = -30. + 0.1 * dindgen(101)

  ;; Model
  case flg of 
     0: begin ;; HRH07
        qso_lf = x_qso_hrh07_lf(z, MODE=-1, M_AB_grid=M_B)
        M1450 = M_B + 0.71 ;; Standard correction (not quite exact)
        qso_lf_eval = interpol(qso_lf, M1450, M1450_EVAL)
     end
     1: begin ;; Willott et al. 2010
        Mstar_1450 = -25.13
        beta = -2.81 ; Bright end
        alpha = -1.5 ; Faint end (assumed)
        Phi_star = 1.14d-8 ; Mpc^-3 mag^-1
        k = -0.47  ; Redshift evolution

        qso_lf_eval = dblarr(n_elements(M1450_EVAL))
        faint = where(M1450_EVAL GE Mstar_1450, nfaint, complement=bright, ncomplement=nbright)
        if nfaint GT 0 then $
           qso_lf_eval[faint] = 10.^(k*(z-6)) * Phi_star / $
                                10.^(0.4 * (alpha+1) * (M1450_EVAL[faint] - Mstar_1450) )
        if nbright GT 0 then $
           qso_lf_eval[bright] = 10.^(k*(z-6)) * Phi_star / $
                                10.^(0.4 * (beta+1) * (M1450_EVAL[bright] - Mstar_1450) )
     end
     else: stop
  endcase

  return, qso_lf_eval

end
