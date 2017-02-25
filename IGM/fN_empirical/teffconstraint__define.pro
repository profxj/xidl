;+ 
; NAME: 
; teffconstraint__define   
;    Version 1.1
;
; PURPOSE:
;    Generate a structure for tau_eff constraints (e.g. D_A) 
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
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
; REVISION HISTORY:
;   July-2011 Written by JXP
;-
;------------------------------------------------------------------------------
pro teffconstraint__define

  tmp = {teffconstraint, $
         ref: '', $              ;; Reference
         cosm: '', $             ;; Cosmology (likely irrelevant)
         type: '', $             ;; e.g. D_A, tau_lyman
         comment: '', $             ;; Fix the type, e.g.  l(X)_tau>2
         z_em: 0., $             ;; Emission Redshift  (sets the Lyman series)
         z_teff: 0., $           ;; Redshift of the evaluation
         NHI_mnx: dblarr(2), $   ;; N_HI to evaluate tau_eff over
         teff: 0., $               ;; tau_eff
         sig_teff: 0. $  ;; sigma(teff)
        }

  return
end
