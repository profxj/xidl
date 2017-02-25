;;
;; Structure holding a MFP constraint
;;
pro betaconstraint__define

  tmp = {betaconstraint, $
         ref: '', $              ;; Reference
         cosm: '', $             ;; Cosmology (e.g. WMAP05)
         type: '\beta', $        ;; Type
         comment: '', $          ;; 
         NHI: fltarr(2), $      ;; Interval where beta is measured
         z_beta: 0., $           ;; Redshift of the evaluation
         beta: 0., $             ;; beta
         sig_beta: 0. $          ;; sigma(beta)
        }

  return
end
