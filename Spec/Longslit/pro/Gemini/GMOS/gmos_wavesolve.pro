PRO GMOS_WAVESOLVE, arcfile, slitfile, wavefile
  
  t0 = systime(1)
  
;----------im
; Read in slit structure 
  ;;slitmask = mrdfits(slitfile, 0) 
  tset_slits = mrdfits(slitfile, 1)
  ystruct = mrdfits(slitfile, 2)
  tslitmask = long_slits2mask(tset_slits) ;; tranformed frame slitmask
  
  nx = tset_slits[0].dims[0]
  ny = tset_slits[0].dims[1]
  idim = size(tset_slits[0].COEFF, /dim)
  nslit = idim[1]
  yarr = findgen(ny) ## replicate(1.0, nx)
  xarr = findgen(nx) # replicate(1.0, ny)


  long_proc, arcfile, arcimg, arcivar, hdr = hdr, /TRANSFORM
  wstruct = long_wstruct(hdr[*, 0], linelist = linelist, reid_file = reid_file)
  qafile = repstr(wavefile, '.fits', '.ps')
  savefile = repstr(wavefile, '.fits', '.sav')
  xfit = long_waveimg(arcimg, arcivar, tset_slits, wstruct, savefile $
                      , fwhmset = fwhmset, qafile = qafile, box_rad = 2.0 $
                      , FWCOEFF = 2)
  ;; This is for the old slits
  ;;wave_slits = dblarr(ny, nslit)
  ;;yvec = dindgen(ny)
  ;;twaveimg = dblarr(nx, ny)
  ;;FOR islit = 0L, nslit-1L DO BEGIN
  ;;   wave_slit_temp = x_calcfit(yvec, fitstr = xfit[islit])
  ;;   ;; convert to vaccuum
  ;;   airtovac, wave_slit_temp
  ;;   wave_slits[*, islit] = wave_slit_temp
  ;;   ipix = WHERE(tslitmask EQ (islit+1L))
  ;;   twave_temp =  x_calcfit(yarr[ipix], fitstr = xfit[islit])
  ;;   ;; convert to vacuum
  ;;   airtovac, twave_temp
  ;;   twaveimg[ipix] = twave_temp
  ;;ENDFOR
  pixset = long_wavepix(arcimg, tset_slits, fwhm = fwhmset.MEDIAN $
                        , box_radius = wstruct.radius $
                        , sig_thresh = wstruct.sig_wpix $
                        , pkwdth = wstruct.pkwdth $
                        , TOLER = wstruct.TOLER, CHK = CHK)
  tpiximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                           , waveimg = twaveimg)
  uwaveimg = gmos_trnimg1to3(twaveimg, tset_slits)
  upiximg = gmos_trnimg1to3(tpiximg,tset_slits)
  
  ;waveimg = waveimg*(slitmask GT 0.0)
  ;;; Now interpolate in the gaps
  ;FOR islit = 0L, nslit-1L DO BEGIN
  ;   gap_pix = WHERE(gap_mask AND slitmask EQ (islit + 1L))
  ;   gd_pix  = WHERE(gap_mask EQ 0 AND slitmask EQ (islit + 1L))
  ;   coeff = poly_fit(ymarr[gd_pix], waveimg[gd_pix], 3)
  ;   waveimg[gap_pix] = poly(ymarr[gap_pix], coeff)
  ;ENDFOR

  mwrfits, uwaveimg, wavefile, /create
  mwrfits, upiximg, wavefile
  mwrfits, fwhmset, wavefile
;;  mwrfits, wave_slits, wavefile
  splog, 'Elapsed time = ', systime(1)-t0, ' sec'
  
RETURN
END
