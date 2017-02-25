pro powerfnstrct__define

  tmp = {powerfnstrct, $
        zmnx: fltarr(2),    $    ;; Range of redshift where the power-laws apply
        npivot: 0, $             ;; Number of pivots
        pivots: dblarr(15), $    ;; log clm where the power-laws pivot
        fn_pivot: fltarr(15), $  ;; Value of log f(N,X) at each pivot
        beta: fltarr(15), $      ;; Power-law for f(N,X)
        zpivot: fltarr(15), $    ;; z pivot in the power-law
        gamma: fltarr(15), $     ;; Power-law exponent for redshift evolution (beware of dX/dz)
        cosmology: fltarr(5)  $  ;; H_0, Omega_m, Omega_L
        }

  return
end
