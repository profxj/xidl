;+ 
; NAME: 
; llsconstraint__define   
;    Version 1.1
;
; PURPOSE:
;    Generate a structure for LLS constraints on f(N,X)
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
pro llsconstraint__define

  tmp = {llsconstraint, $
         ref: '', $              ;; Reference
         cosm: '', $             ;; Cosmology (e.g. WMAP05)
         type: '', $             ;; Fix the type, e.g.  l(X)_tau>2
         comment: '', $             ;; Fix the type, e.g.  l(X)_tau>2
         z_lls: 0., $            ;; Redshift of the evaluation
         tau_lim: 0., $          ;; Minimum tau
         lX: 0., $                ;; l(X)
         sig_lX: 0.$               ;; sigma(lX)
        }

  return
end
