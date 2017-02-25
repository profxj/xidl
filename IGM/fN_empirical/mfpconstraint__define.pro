;;
;; Structure holding a MFP constraint
;;
pro mfpconstraint__define

  tmp = {mfpconstraint, $
         ref: '', $              ;; Reference
         cosm: '', $             ;; Cosmology (e.g. WMAP05)
         type: '\lmfp', $             ;; Fix the type, e.g.  l(X)_tau>2
         comment: '', $             ;; Fix the type, e.g.  l(X)_tau>2
         z_mfp: 0., $            ;; Redshift of the evaluation
         mfp: 0., $               ;; MFP
         sig_mfp: 0. $  ;; sigma(MFP)
        }

  return
end
