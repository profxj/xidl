

;; H+K setup
wvmnx = [1.3d, 2.1d] ;; expand the ragne to get slightly more coverage
;;dlam = (2.30-2.0)/1024.0d
lam_cen = 1.6490627
;;lam_cen = 2.0 + (2.30-2.0)/2.0d
;; According to manual, assume resolution is actually 2 pixels
dlam = 3.4778201/1d4
;; arc lines are about 4.5 pixels
resolution = (lam_cen/(4.5*dlam))
;;2200.0 ;; lam_cen/(2.5*dlam)
flgd = 0
linefile = 'SOFI_modelsky_OH_linelist_H.lst'
outfile = 'SOFI_modelsky_OH_linelist_H.fits'
plate_scale = 0.288d
fnslit = 0.6/plate_scale ;; slits are 2x plate scale
;;fnslit = 3.0d
;;pkwdth = 3.0 ;;1.3*fnslit
;;toler = fnslit/2.0d > 2.0d
;;toler = 3.0
pkwdth = 1.3*plate_scale
toler = fnslit/2.0d > 2.0d
ICLSE = 2.0
thin = 1
fweight = 0
nsig = 4.0d
;;T_BB = 273.0
;;T_BB = 250.0
;;wave_water = 2.33
nearir_modelsky_linelist, resolution, linefile, outfile $
                          , wvmnx = wvmnx, dlam = dlam, flgd = flgd $
                          , pkwdth = pkwdth, toler = toler, thin = thin $
                          , fweight = fweight, nsig = nsig, T_BB = T_BB $
                          , wave_water = wave_water, ICLSE = ICLSE
END
