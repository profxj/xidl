pro splinefnstrct__define

  tmp = {splinefnstrct, $
        zmnx: fltarr(2),    $    ;; Range of redshift where the spline applies
        npivot: 0, $             ;; Number of spline points
        pivots: dblarr(20), $    ;; log clm where the power-laws pivot
        fn_pivot: dblarr(20), $  ;; Value of log f(N,X) at each pivot
        splin: dblarr(20), $  ;; Value of spline initialization
        zpivot: fltarr(20), $    ;; z pivot in the power-law
        gamma: fltarr(20), $     ;; Power-law exponent for redshift evolution (beware of dX/dz)
        cosmology: fltarr(5)  $  ;; H_0, Omega_m, Omega_L
        }

  return
end
