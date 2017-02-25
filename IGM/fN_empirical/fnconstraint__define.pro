;;
;; Structure holding a specific set of f(N) constraints from the
;; literature
;;
pro fnconstraint__define

  tmp = {fnconstraint, $
         ref: '', $              ;; Reference
         cosm: '', $             ;; Cosmology (e.g. WMAP05)
         type: '', $             ;; e.g. Lya Forest, SLLS, DLA
         comment: '', $          ;; Comment
         zeval: 0., $            ;; Redshift of the evaluation
         npt: 0, $               ;; Number of values 
         NHI: fltarr(100), $     ;; log N_HI values
         bins: fltarr(100,2), $  ;; Bin points in N_HI
         DX: 0., $               ;; Total Delta X in the analysis
         DN: fltarr(100), $      ;; Delta N for each f(N) bin
         fN: fltarr(100), $      ;; log f(N_HI,X)
         sig_fN: fltarr(100,2) $ ;; sigma[log f(N_HI,X)]  (hi/low)
        }

  return
end
