PRO GMOS_NS_EXTRACT, scifile, darkfile, slitfile, wavefile, outfile $
                     , EDG_PIX = EDG_PIX1, SIGREJ = SIGREJ1 $
                     , NOVACHELIO = NOVACHELIO $
                     , NOZAP = NOZAP, NOTRAPZAP = NOTRAPZAP
  
  ;;IF NOT KEYWORD_SET(BOX_RAD1) THEN box_rad = 2.5D $
  ;;ELSE BOX_RAD = BOX_RAD1
  IF KEYWORD_SET(SIGREJ1) THEN SIGREJ = SIGRE1 $
  ELSE sigrej = 5.0
  
  ;; proc image with dark and bias subtraction
  long_proc, scifile, sciimg, sciivar, hdr = hdr, biasfile = darkfile
  grating = strcompress(sxpar(hdr[*, 0], 'GRATING'), /rem)
  IF KEYWORD_SET(EDG_PIX) THEN EDG_PIX = EDG_PIX1 $
  ELSE BEGIN 
     IF strmatch(grating, '*R150*') THEN EDG_PIX = 1.0 $
     ELSE EDG_PIX = 1.0
  ENDELSE
  slitmask = mrdfits(slitfile, 0)
  tset_slits = mrdfits(slitfile, 1)
  ;; Create some images we will need later. These are in transformed frame
  tximg = long_slits2x(tset_slits,TOL_EDG=EDG_PIX,EDGMASK=TEDGMASK)
  ;; Map to untransformed frame
  ;;ximg = gmos_trnimg1to3(tximg, tset_slits)
  edgmask = gmos_trnimg1to3(tedgmask, tset_slits)

  ystruct = mrdfits(slitfile, 2)
  nx = tset_slits[0].DIMS[0]
  ny = tset_slits[0].DIMS[1]
  nxby3 = nx/3
  idim = size(tset_slits[0].COEFF, /dim)
  nslit = idim[1]
  slit_width = tset_slits[1].coeff[0, *] - tset_slits[0].coeff[0, *]
  med_width = djs_median(tset_slits[1].coeff[0, *] - tset_slits[0].coeff[0, *])
  ;; Determine the binning
  gmos_binning,nx,ny,specbin=specbin,spatbin=spatbin,BTAG=BTAG
  IF strmatch(grating, '*R150*') THEN BEGIN
     ;; For the R150 grating one resolution element is 5 pixels
     ;; so we bin in the extraction by a factor of 3
     fnres = 2.0D
     nopt =  long(2400.0d/fnres)
     cdelt_res = 3.515d
     cdelt = cdelt_res*double(fnres)
     wave_opt = 3000.0D + dindgen(nopt)*cdelt
  ENDIF ELSE BEGIN
     ;; For the R400 grating one resolution element is 5 pixels
     ;; so we bin in the extraction by a factor of 3
     fnres = 2.0D 
     nopt = long(4500.0d/fnres)
     cdelt_res = 1.3710D
     cdelt = cdelt_res*double(fnres)
     wave_opt = 4000.0D + dindgen(nopt)*cdelt
  ENDELSE
  proto = create_struct('WAVE_BOX', dblarr(ny) $     ; boxcar wavelengths
                        , 'FLUX_BOX', fltarr(ny) $   ; boxcar flux
                        , 'IVAR_BOX', fltarr(ny) $   ; boxcar inverse var
                        , 'SKY_BOX', fltarr(ny) $    ; boxcar sky
                        , 'RN2_BOX', fltarr(ny) $    ; boxcar RN^2
                        , 'PIX_BOX', fltarr(ny) $    ; total pixels in boxcar
                        , 'MASK_BOX', bytarr(ny) $   ; boxcar mask
                        , 'WAVE_OPT', wave_opt $
                        , 'FLUX_OPT',fltarr(nopt) $
                        , 'IVAR_OPT',fltarr(nopt) $
                        , 'SKY_OPT',fltarr(nopt)  $
                        , 'RN2_OPT',fltarr(nopt)  $
                        , 'PIX_OPT_TOT', fltarr(nopt)  $
                        , 'PIX_OPT_USE', fltarr(nopt)  $
                        , 'MASK_OPT',bytarr(nopt) $
                        , 'XPOS', fltarr(ny) $   ; midpoint of trace x
                        , 'YPOS', fltarr(ny) $   ; midpoint of trace y
                        , 'SEDG_L', fltarr(ny) $ ; slit edges in untrans frm
                        , 'SEDG_R', fltarr(ny) $  
                        , 'BOX_RAD', 0.0 $ ;; BOX_RAD $
                        , 'EDG_PIX', EDG_PIX $  
                        , 'EXPTIME', 0.0) ; boxcar radius
;; left_ap is SEDG_L + EDG_PIX, righ_ap is SEDG_R - EDG_PIX
  objstruct = replicate(proto, nslit)
  objstruct = struct_addtags(objstruct, ystruct)
  ;; read in the wavelength map
  waveimg    = mrdfits(wavefile, 0)
  piximg     = mrdfits(wavefile, 1)
  fwhmset    = mrdfits(wavefile, 2)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;  vacuum and heliocentric correction
  IF NOT keyword_set(NOVACHELIO) THEN BEGIN
     iwv=WHERE(waveimg GT 0)
     tmp_waveimg=waveimg[iwv]
     airtovac, tmp_waveimg
     waveimg[iwv]=tmp_waveimg
     if n_elements(size(hdr,/dime)) EQ 1 then idx = 0 else idx = 1
     mjd = double(sxpar(hdr[*, idx], 'MJD-OBS'))  + 2400000.5D
     equinox = double(sxpar(hdr[*, idx], 'EQUINOX'))
     radeg = double(sxpar(hdr[*, 0], 'RA'))
     decdeg = double(sxpar(hdr[*, 0], 'DEC'))
     helio = (-1.0D)*x_keckhelio(radeg, decdeg, equinox, jd = mjd)
     print, 'gmos_extract: heliocentric correction :', $
            helio, ' km/s', format = '(a,f8.3,a)'
     hel_corr = sqrt( (1.d + helio/299792.458d) / $
                      (1.d - helio/299792.458d))
     waveimg = waveimg*hel_corr
  ENDIF
  IF KEYWORD_SET(darkfile) THEN darkmask = mrdfits(darkfile, 1) $
  ELSE message, 'ERROR: Missing darks'
  badpix = WHERE(darkmask EQ 0, nbad)
  IF nbad GT 0 THEN sciivar[badpix] = 0.0
;; ANOD =  Actual number of non-shuffled exposures begun  
;; BNOD =  Actual number of shuffled exposures begun    
  ANOD = long(sxpar(hdr[*, 0], 'ANODCNT'))
  BNOD = long(sxpar(hdr[*, 0], 'BNODCNT'))
  exptime = float(sxpar(hdr[*, 0], 'EXPOSURE'))
  objstruct.EXPTIME = exptime
;; A region is nbxy3:2*nxby3-1
;; B region is 0:nxby3
  nod_scale = double(ANOD)/double(BNOD)
  sciimg_A  = sciimg[nxby3:2*nxby3-1L, *]
  sciivar_A = sciivar[nxby3:2*nxby3-1L, *]
  sciimg_B  = float(nod_scale)*sciimg[0:nxby3-1L, *]
  sciivar_B = sciivar[0:nxby3-1L, *]/nod_scale^2
  img_minsky  = sciimg_A - sciimg_B
  mask = (sciivar_A GT 0.0 AND sciivar_B GT 0.0)
  var_minsky = float(mask)*(1.0/(sciivar_A + (sciivar_A EQ 0.0)) + $
                            1.0/(sciivar_B + (sciivar_B EQ 0.0)))
  ivar_minsky = float(mask)/(var_minsky + (var_minsky EQ 0.0))
  slitmask = slitmask[nxby3:2*nxby3-1L, *]
  waveimg  = waveimg[nxby3:2*nxby3-1L, *]
  edgmask = edgmask[nxby3:2*nxby3-1L, *]
  ;; Cosmic ray reject
  IF NOT KEYWORD_SET(NOZAP) THEN BEGIN
     ;;   Take the PSF width to be that of the spectral direction.  
     ;;   This prevents the routine from rejecting sky lines
     ;;   This is a description of the 3x3 core of the 2D PSF for reject_cr.pro
     ;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1] 
     ;;    PSFVALS[0]          1.   PSFVALS[0]
     ;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
     fwhm_vec = [djs_median(fwhmset.median)] ;;, med_width]
     fwhm = max(fwhm_vec)
     fudge = 0.5d
     sigma_psf = fwhm/2.35482D*fudge
     psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
     ncr = 0
     ;; interpolate in y-dir since most features are in x
     ivar_interp = djs_maskinterp(ivar_minsky, mask EQ 0, iaxis = 1)
     ivar_smth = smooth(ivar_interp, [1, fwhm], /edge_truncate)
     ibad = WHERE(mask EQ 0, nbad)
     IF nbad GT 0 THEN ivar_interp[ibad] = ivar_smth[ibad]
     ivar_smth = smooth(ivar_interp, [1, fwhm], /edge_truncate)
     ;;ivar_smth = median(ivar_smth, fwhm)*(ivar_minsky GT 0.0)
     djs_iterstat, img_minsky[where(mask AND slitmask)], sigma = sigma
     crmaskA  = psf_reject_cr(img_minsky > (-3*sigma), ivar_smth, psfvals $
                              , satmask = (img_minsky GT 1d5 OR mask EQ 0));; $
     ;;OR slitmask EQ 0))
     ;; update mask with img A crs
     mask = mask*(crmaskA EQ 0)
     ivar_interp = djs_maskinterp(ivar_minsky, mask EQ 0, iaxis = 1)
     ivar_smth1 = smooth(ivar_interp, [1, fwhm], /edge_truncate)
     ibad = WHERE(mask EQ 0, nbad)
     IF nbad GT 0 THEN ivar_interp[ibad] = ivar_smth1[ibad]
     ivar_smth = smooth(ivar_interp, [1, fwhm], /edge_truncate)
     crmaskB = psf_reject_cr(-img_minsky > (-3*sigma), ivar_smth, psfvals $
                             , satmask = (img_minsky GT 1d5 OR mask EQ 0)) ;; $
     ;;OR slitmask EQ 0))
     ;; update mask with img B crs
     mask = mask*(crmaskB EQ 0)
     delvarx, ivar_smth, ivar_interp
  ENDIF
  IF NOT KEYWORD_SET(NOTRAPZAP) THEN BEGIN
     ;; Now filter our residual charge horizontal traps that do not
     ;; not get removed in dark masking or CR zapping
     smash = djs_avsigclip(img_minsky*mask, 1, inmask = (mask EQ 0) $
                           , outmask = outmask)
     sig_res_med  =  round(ny/300L)
     sig_res_spec = 2*sig_res_med
     smash_med = djs_median(smash, width = sig_res_med, boundary = 'reflect')
     nhalf =  long(sig_res_spec)*4L
     xkern = dindgen(2*nhalf+1)-nhalf
     kernel = gauss1(xkern, [0.0, sig_res_spec, 1.0])
     smash_smth = convol(smash_med, kernel, /edge_truncate)
     djs_iterstat, smash-smash_smth, sigma = sig_sm
     bad1 = where(abs(smash-smash_smth) GT 7.0*sig_sm, nbad1)
     dummask1 = lonarr(ny) + 1L
     IF nbad1 GT 0 THEN dummask1[bad1] = 0L
     ;; Grow this mask by 1 pixel
     dummask = long_grow_mask(dummask1, 1)
     bad = where(dummask EQ 0, nbad)
     ;;ymask = 0*mask + 1L
     IF nbad GT 0 THEN mask[*, bad] = 0
     ivar_minsky = ivar_minsky*mask
  ENDIF
  
  ;; Expand slits in transformed frame
  traceset2xy, tset_slits[0], y_tr, xl_tr
  traceset2xy, tset_slits[1], y_tr, xr_tr
  ;; Map transformed frame slits to untransformed frame
  gmos_trnxy1to3, reform(xl_tr, ny*nslit), reform(y_tr, ny*nslit) $
                  , x3 = xl1, y3 = yrow1,BTAG=BTAG 
  gmos_trnxy1to3, reform(xr_tr, ny*nslit), reform(y_tr, ny*nslit) $
                  , x3 = xr1, y3 = yrow1,BTAG=BTAG
;; now transform everything over to 0-768 frame
  xl1  = reform(xl1 - float(nxby3), ny, nslit)
  xr1  = reform(xr1 - float(nxby3), ny, nslit)
  yrow1= reform(yrow1,ny,nslit)
  xl = fltarr(ny,nslit)
  xr = fltarr(ny,nslit)
  ;; Now interpolate the slit boundaries onto a fixed grid in y
  yvec=findgen(ny)
  FOR islit=0L,nslit-1L DO BEGIN
     ;; make sure this is monotonic
     isort=sort(yrow1[*,islit])
     xl_temp = interpol(xl1[isort,islit],yrow1[isort,islit],yvec)
     xl[isort,islit] = xl_temp 
     xr_temp = interpol(xr1[isort,islit],yrow1[isort,islit],yvec)
     xr[isort,islit] = xr_temp 
  ENDFOR
  yrow=y_tr ;; yrow is now a fixed grid (same in both frames)
  ;; trace position is midpoint of slit, extract asymmetric boxcars
  ;; EDG_PIX pixels away from slit edge on each side
  trace = (xl + xr)/2.0D
  ap_L = xl + EDG_PIX
  ap_R = xr - EDG_PIX

  skysub_mask = fltarr(nxby3, ny)
  sky_img =  fltarr(nxby3, ny)
  sign = fltarr(nslit)
  FOR islit = 0L, nslit-1L DO BEGIN
     ipix = WHERE(slitmask EQ (islit + 1L), nonmask)
     IF nonmask EQ 0 THEN BEGIN
        splog, 'WARNING: slit #', islit + 1L, ' was not found'
        splog, 'It could be outside of the NS region or there is an error'
        sign[islit] = 0.0
     ENDIF ELSE BEGIN
        sign[islit] = (objstruct[islit].SKYFLAG EQ 1) ? 1.0 : -1.0
        skysub_mask[ipix] = sign[islit]
     ENDELSE
 ENDFOR
  rn_img =  long_rdnoiseimg(nx, ny, hdr)
  rn2_img   = rn_img[nxby3:2*nxby3-1L, *]^2 +  $
              float(nod_scale)^2*rn_img[0:nxby3-1L, *]^2
  iskya =  WHERE(skysub_mask EQ 1)
  iskyb =  WHERE(skysub_mask EQ -1)
  sky_img[iskya]  = sciimg_A[iskya]
  sky_img[iskyb]  = sciimg_B[iskyb]
  
  sign_arr = replicate(1.0, ny) # sign 
  newmask = (ivar_minsky GT 0.0) ;;AND (slitmask GT 0.0)
  rawflux = extract_asymbox2(img_minsky*newmask,ap_L,ap_R)
  rawsky  = extract_asymbox2(sky_img*newmask, ap_L,ap_R)
  rawrn2  = extract_asymbox2(rn2_img*newmask, ap_L,ap_R)
  box_denom = extract_asymbox2((waveimg GT 0.0), ap_L,ap_R)
  wave_now =  extract_asymbox2(waveimg, ap_L,ap_R)/ $
              (box_denom + (box_denom EQ 0))
  ;; interpolate over places where wavelengths are zero
  FOR islit=0L,nslit-1L DO BEGIN
     ibad = WHERE(wave_now[*,islit] LE 0.0,nwvbad $
                  ,COMPLEMENT=igood,NCOMP=nwvgood)
     IF nwvbad GT 0 THEN BEGIN
        isort = sort(yvec[igood])
        wave_now[ibad,islit]= $
           interpol(wave_now[igood[isort],islit],yvec[igood[isort]],yvec[ibad])
     ENDIF
  ENDFOR
  objstruct.WAVE_BOX = wave_now
  ;;objstruct.WAVE = wave_slits
  varimg = 1.0D/(ivar_minsky + (ivar_minsky EQ 0))
  var_box = extract_asymbox2(varimg, ap_L,ap_R)
  ivar_box = 1.0/(var_box + (var_box EQ 0))
  pixtot  = extract_asymbox2(float(newmask*0 + 1.0D), ap_L,ap_R)
  pixgood = extract_asymbox2(float(newmask), ap_L,ap_R)
  ;;scale_ratio = pixtot/(pixgood + (pixgood EQ 0.0))
  ;; What is the proper boxcar?? I think it is scale_ratio =1
  scale_ratio = 1.0d
  objstruct.MASK_BOX = $
     extract_asymbox2(float(newmask EQ 0), ap_L,ap_R) NE pixtot
  objstruct.FLUX_BOX = rawflux*scale_ratio*objstruct.MASK_BOX*sign_arr
  objstruct.IVAR_BOX = ivar_box*objstruct.MASK_BOX/scale_ratio^2
  objstruct.SKY_BOX  = scale_ratio*rawsky
  objstruct.RN2_BOX  = scale_ratio*rawrn2
  objstruct.PIX_BOX = pixgood
  objstruct.SEDG_L = xl
  objstruct.SEDG_R = xr
  objstruct.XPOS = trace
  objstruct.YPOS = yrow
  
  xarr=lindgen(nxby3) # replicate(1L,ny)
  ;;optmask = lonarr(nxby3, ny) + 1
  ;; Optimal extraction 
  FOR islit=0L,nslit-1L DO BEGIN
     slitpix = WHERE(slitmask EQ (islit+1L),nslp)
     IF nslp EQ 0 THEN CONTINUE
     ;; Current min amd max wavelength for this slit
     mwv=minmax(waveimg[slitpix])
     min_wave=mwv[0] - cdelt
     max_wave=mwv[1] + cdelt
     dwave = max_wave - min_wave
     mmx = minmax(xarr[slitpix])
     min_x=mmx[0]
     max_x=mmx[1]
     ;;  slit counter
     IF NOT KEYWORD_SET(SILENT) THEN $
        print, format = '("EXTRACTING SLIT ",i2," of ",i2,a1,$)' $
               , islit+1L, nslit, string(23b)
     FOR iwv=0L,nopt-1L DO BEGIN
        IF wave_opt[iwv] LT min_wave OR wave_opt[iwv] GT max_wave THEN CONTINUE
        min_y = floor((wave_opt[iwv] - 50.0*cdelt_res - min_wave)/dwave*ny) > 0
        max_y = ceil((wave_opt[iwv]  + 50.0*cdelt_res - min_wave)/dwave*ny) $
                < (ny-1L)
        nc = max_x - min_x + 1L
        nr = max_y - min_y + 1L
        ipix = lindgen(nc) # replicate(1, nr) + $
               replicate(1, nc) # (min_y + lindgen(nr))*nxby3 + min_x
        ibin=WHERE(slitmask[ipix] EQ (islit+1L) AND $
                   edgmask[ipix] EQ 0 AND $
                   waveimg[ipix] GE (wave_opt[iwv] - cdelt/2.0D) AND $
                   waveimg[ipix] LT (wave_opt[iwv] + cdelt/2.0D) AND $
                   newmask[ipix] GT 0,nbin)
        objstruct[islit].WAVE_OPT = wave_opt
        IF nbin GT 0 THEN BEGIN
           djs_iterstat, img_minsky[ipix[ibin]], sigrej = 3.0, mean = flux_opt $
                         , mask = mask_rej
           objstruct[islit].FLUX_OPT[iwv] = flux_opt*sign[islit]
           pix_opt=total(mask_rej)
           objstruct[islit].PIX_OPT_USE[iwv]=pix_opt
           objstruct[islit].PIX_OPT_TOT[iwv] = float(nbin)
           var_opt = total(varimg[ipix[ibin]]*float(mask_rej))/pix_opt^2
           objstruct[islit].IVAR_OPT[iwv] = $
              (pix_opt GT 0)/(var_opt + (var_opt EQ 0))
           objstruct[islit].SKY_OPT[iwv] = $
              total(sky_img[ipix[ibin]]*float(mask_rej))/pix_opt
           objstruct[islit].RN2_OPT[iwv] = $
              total(rn2_img[ipix[ibin]]*float(mask_rej))/pix_opt
           objstruct[islit].MASK_OPT[iwv] = (pix_opt GT 0)
           ;; mask these rejected pixels now in the mask
           mask[ipix[ibin]] = mask_rej
        ENDIF ELSE CONTINUE
     ENDFOR
     ;; Now compute the average pix_opt_tot for each aperture and scale 
     ;; everything appropriately to match total boxcar flux
     djs_iterstat, objstruct[islit].PIX_OPT_TOT, sigrej = 3.0 $
                   , mean = pix_opt_tot
     objstruct[islit].FLUX_OPT =  objstruct[islit].FLUX_OPT*pix_opt_tot
     objstruct[islit].IVAR_OPT =  objstruct[islit].IVAR_OPT/pix_opt_tot^2
     objstruct[islit].SKY_OPT  =  objstruct[islit].SKY_OPT*pix_opt_tot
     objstruct[islit].RN2_OPT  = objstruct[islit].RN2_OPT*pix_opt_tot
  ENDFOR
  fxhmake, hdr1, float(img_minsky), /EXTEND, /INIT
  sxdelpar, hdr1, 'END'
  hdr_old = hdr[*, 0]
  sxdelpar, hdr_old, 'NAXIS'
  sxdelpar, hdr_old, 'NAXIS2'
  sxdelpar, hdr_old, 'BZERO'
  sxdelpar, hdr_old, 'BSCALE'
  sxdelpar, hdr_old, 'SIMPLE'
  sxdelpar, hdr_old, 'BITPIX'
  sxdelpar, hdr_old, 'EXTEND'
  hdr_out = [hdr1, hdr_old]
  ;; Fix headers in case there are things we need in o
  mwrfits, float(img_minsky), outfile, hdr_out, /create
  mwrfits, float(ivar_minsky*mask), outfile
  mwrfits, float(skysub_mask), outfile
  mwrfits, double(waveimg),outfile
  mwrfits, objstruct, outfile
  mwrfits, float(sciimg), outfile
  mwrfits, float(sciivar), outfile
  RETURN
END
