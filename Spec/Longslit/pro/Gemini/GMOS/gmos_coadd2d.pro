PRO GMOS_COADD2D, infiles, slitfiles, slitno, outfile2d, lambda_ref $
                  , SIGREJ = SIGREJ, outfile1d = outfile1d
  
  head0 = xheadfits(infiles[0])
  IF NOT KEYWORD_SET(SIGREJ1) THEN sigrej = 3.0 $
  ELSE SIGREJ = SIGREJ1
  nap = n_elements(slitno)
  nfiles = n_elements(infiles)
  nimgs = nap*nfiles
  xsize = lonarr(nimgs)
  lambda_ref1 = 7180.0d
;; First pass read in object structures. 
;; determine how large an array we need
  FOR ifile = 0L, nfiles-1L DO BEGIN
     IF ifile EQ 0 THEN BEGIN
        tset_slits = mrdfits(slitfiles[ifile], 1, /silent)
        nx = tset_slits[0].DIMS[0]
        ny = tset_slits[0].DIMS[1]
        nxby3 = nx/3
        idim = size(tset_slits[0].COEFF, /dim)
        nslit = idim[1]
        xpix = findgen(nxby3) # replicate(1.0, ny)
        ypix = replicate(1.0, nxby3) # findgen(ny)
     ENDIF
     FOR iap = 0L, nap-1L DO Begin
        obj1 = mrdfits(infiles[ifile], 4, hdr, /silent)
        inow = WHERE(obj1.SLIT EQ slitno[iap], nthis)
        IF nthis EQ 0 THEN message, 'Problem with input slitno. Could not find' 
        IF NOT KEYWORD_SET(OBJSTRUCT) THEN objstruct = obj1[inow] $
        ELSE objstruct = [objstruct, obj1[inow]]
     ENDFOR
  ENDFOR
  exptime = total(objstruct.EXPTIME)
  iref = 0L
  objstruct_new = objstruct
  slit_extent = objstruct.SEDG_R - objstruct.SEDG_L + 1L
  nxcut1 = max(slit_extent)
  min_x = floor(min(objstruct[iref].SEDG_L - nxcut1/2L))
  max_x = ceil(max(objstruct[iref].SEDG_R + nxcut1/2L))
  nxcut = max_x - min_x + 1L
  ;; size must be even
;IF nxcut MOD 2 EQ 0 THEN BEGIN
;   IF min_x GT 1 THEN BEGIN
;      nxhalf = nxcut/2L
;      nxcut  = nxcut+1L
;      min_x = min_x - 1L ;; extend one in min to make size odd
;   ENDIF ELSE BEGIN
;      nxhalf = nxcut/2L
;      nxcut  = nxcut+1L
;      max_x = max_x + 1L 
;   ENDELSE
;ENDIF ELSE nxhalf = (nxcut - 1L)/2L
  img_arr  = fltarr(nxcut, ny, nimgs)
  var_arr = fltarr(nxcut, ny, nimgs)
  waveimg_arr = fltarr(nxcut, ny, nimgs)
  skysub_mask_arr = fltarr(nxcut, ny, nimgs)
  slitmask_arr = fltarr(nxcut, ny, nimgs)
;; What is the reference wavelength 
  img_minsky = mrdfits(infiles[iref], 0, hdr, /silent)
  slitmask = mrdfits(slitfiles[iref], 0)
  slitmask = slitmask[nxby3:2*nxby3-1L, *]
  ivar        = mrdfits(infiles[iref], 1, /silent)
  skysub_mask = mrdfits(infiles[iref], 2, /silent)
  waveimg     = mrdfits(infiles[iref], 3, /silent)
  
  min_diff = min(abs(objstruct[iref].WAVE_BOX - lambda_ref1), pix_y_ref)
  lambda_ref = objstruct[iref].WAVE_BOX[pix_y_ref]
  pix_x_ref = interpol(objstruct[iref].XPOS, objstruct[iref].YPOS, pix_y_ref)
  img_arr[*, *, iref] = img_minsky[min_x:max_x, *]
  slitmask_arr[*, *, iref] = slitmask[min_x:max_x, *] EQ objstruct[iref].SLIT
  var = (ivar GT 0)/(ivar + (ivar LE 0))
  var_arr[*, *, iref] = var[min_x:max_x, *]
  skysub_mask_arr[*, *, iref] = skysub_mask[min_x:max_x, *]
  waveimg_arr[*, *, iref] =  waveimg[min_x:max_x, *]
  
  objstruct_new[iref].XPOS = objstruct[iref].xpos - min_x
  objstruct_new[iref].SEDG_L = objstruct_new[iref].SEDG_L - min_x
  objstruct_new[iref].SEDG_R = objstruct_new[iref].SEDG_R - min_x
  FOR ifile = 0L, nfiles-1L DO BEGIN
     img_minsky = mrdfits(infiles[ifile], 0, hdr, /silent)
     slitmask = mrdfits(slitfiles[ifile], 0)
     slitmask = slitmask[nxby3:2*nxby3-1L, *]
     ivar        = mrdfits(infiles[ifile], 1, /silent)
     var =  (ivar GT 0)/(ivar + (ivar LE 0))
     skysub_mask = mrdfits(infiles[ifile], 2, /silent)
     waveimg     = mrdfits(infiles[ifile], 3, /silent)
     FOR iap = 0L, nap-1L DO BEGIN
        objind = ifile*nap + iap
        IF objind EQ iref THEN CONTINUE
        pix_y_shf = round(interpol(objstruct[objind].YPOS $
                                   , objstruct[objind].WAVE_BOX, lambda_ref))
        yshift = pix_y_ref - pix_y_shf 
        pix_x_shf = interpol(objstruct[objind].XPOS, objstruct[objind].YPOS $
                             , pix_y_ref)
        xshift = round(pix_x_ref  - pix_x_shf)
        ;; img
        img_tmp = shift(img_minsky, xshift, yshift)
        IF yshift GT 0 THEN img_tmp[*, 0:(yshift-1L)] = 0 $
        ELSE IF yshift LT 0 THEN img_tmp[*, (ny+yshift):*] = 0
        img_arr[*, *, objind] = img_tmp[min_x:max_x, *]
        ;; slitmask
        slitmask_tmp = shift((slitmask EQ objstruct[objind].SLIT) $
                             , xshift, yshift)
        IF yshift GT 0 THEN slitmask_tmp[*, 0:(yshift-1L)] = 0 $
        ELSE IF yshift LT 0 THEN slitmask_tmp[*, (ny+yshift):*] = 0
        slitmask_arr[*, *, objind] = slitmask_tmp[min_x:max_x, *]
        ;; var
        var_tmp = shift(var, xshift, yshift)
        IF yshift GT 0 THEN var_tmp[*, 0:(yshift-1L)] = 0 $
        ELSE IF yshift LT 0 THEN var_tmp[*, (ny+yshift):*] = 0
        var_arr[*, *, objind] = var_tmp[min_x:max_x, *]
        ;; skysub_mask
        skysub_mask_tmp = shift(skysub_mask, xshift, yshift)
        IF yshift GT 0 THEN skysub_mask_tmp[*, 0:(yshift-1L)] = 0 $
        ELSE IF yshift LT 0 THEN skysub_mask_tmp[*, (ny+yshift):*] = 0
        skysub_mask_arr[*, *, objind] = skysub_mask_tmp[min_x:max_x, *]
        ;; waveimg
        waveimg_tmp = shift(waveimg, xshift, yshift)
        IF yshift GT 0 THEN waveimg_tmp[*, 0:(yshift-1L)] = 0 $
        ELSE IF yshift LT 0 THEN waveimg_tmp[*, (ny+yshift):*] = 0
        waveimg_arr[*, *, objind] = waveimg_tmp[min_x:max_x, *]
        
        objstruct_new[objind].XPOS = objstruct[objind].xpos - min_x + xshift
        objstruct_new[objind].SEDG_L = objstruct_new[objind].SEDG_L - min_x $
                                       + xshift
        objstruct_new[objind].SEDG_R = objstruct_new[objind].SEDG_R - min_x $
                                       + xshift
     ENDFOR
  ENDFOR
  mask_arr = (var_arr GT 0)*(slitmask_arr GT 0) EQ 0
;; Average the images with rejection. Mask convention: 0 = good, 1=bad
  img_avs = djs_avsigclip(img_arr*skysub_mask_arr, 3, sigrej = sigrej $
                          , maxiter = maxiter, inmask = mask_arr $
                          , outmask = outmask_arr)
  weights = img_arr*0 + 1.0 ;; unit weights for now
; Combine the variances
  nused = total(outmask_arr EQ 0, 3)
  varfinal_avs = total(var_arr*(outmask_arr EQ 0), 3)/(nused^2 + (nused EQ 0))
  maskfinal = (total(outmask_arr, 3) NE nimgs)
  ;; Optimally combine the images
  weights = weights*float(outmask_arr EQ 0)
  wght_sum = total(weights, 3)
  imgfinal = total(weights*img_arr*skysub_mask_arr $
                   , 3)/(wght_sum + (wght_sum EQ 0.0))
  wavefinal = total(weights*waveimg_arr, 3)/(wght_sum + (wght_sum EQ 0.0))
  varfinal = total(weights^2*var_arr, 3)/(wght_sum + (wght_sum EQ 0.0))^2

  IF KEYWORD_SET(OUTFILE1D) THEN BEGIN
     ;; Find the two closest aperture boundaries
     ap_L = objstruct_new[0].SEDG_L
     ap_R = objstruct_new[0].SEDG_R
     FOR iobj = 0L, nimgs-1L DO BEGIN
        ap_L = ap_L > objstruct_new[iobj].SEDG_L
        ap_R = ap_R < objstruct_new[iobj].SEDG_R
        ;;objstruct_new[iobj].SEDG_L
     ENDFOR
     yvec = findgen(ny)
     ;; Now also do a 1d boxcar extraction
     rawflux_box = extract_asymbox2(imgfinal*maskfinal, ap_L, ap_R)
     box_denom = extract_asymbox2((wavefinal GT 0.0), ap_L, ap_R)
     wave_box =  extract_asymbox2(wavefinal, ap_L, ap_R)/ $
                 (box_denom + (box_denom EQ 0))
     ;; interpolate over places where wavelengths are zero
     ibad = WHERE(wave_box LE 0.0, nwvbad $
                  , COMPLEMENT = igood, NCOMP = nwvgood)
     IF nwvbad GT 0 THEN BEGIN
        isort = sort(yvec[igood])
        wave_box[ibad] = $
           interpol(wave_box[igood[isort]], yvec[igood[isort]], yvec[ibad])
     ENDIF
     rawvar_box = extract_asymbox2(varfinal, ap_L, ap_R)
     pixtot  = extract_asymbox2(float(maskfinal*0 + 1.0D), ap_L, ap_R)
     pixgood = extract_asymbox2(float(maskfinal), ap_L, ap_R)
     ;; What is the proper boxcar?? I think it is scale_ratio =1
     MASK_BOX = extract_asymbox2(float(maskfinal EQ 0), ap_L, ap_R) NE pixtot
     flux_box =   rawflux_box*MASK_BOX
     var_box  =  rawvar_box*mask_box
     hdr_tmp = head0
     sxaddpar, hdr_tmp, 'NEXP', nimgs
     sxaddpar, hdr_tmp, 'EXPTIME_TOT', exptime
     sxaddpar, hdr_tmp, 'BITPIX', -32
     sxaddpar, hdr_tmp, 'NAXIS', 1
     sxaddpar, hdr_tmp, 'NAXIS1', n_elements(flux_box)
     sxdelpar, hdr_tmp, 'NAXIS2'
     sxdelpar, hdr_tmp, 'BZERO'
     sxdelpar, hdr_tmp, 'BSCALE'
     mwrfits, flux_box, outfile1d, hdr_tmp, /create
     sig = sqrt(var_box)
     mwrfits, sig, outfile1d
     mwrfits, wave_box, outfile1d
     print, 'gmos_coadd2d: Final file 1-d file is ', outfile1d
  ENDIF
  ;;write out to file
  fxhmake, hdr1, float(imgfinal), /EXTEND, /INIT
  sxdelpar, hdr1, 'END'
  hdr_old = head0[*, 0]
  sxdelpar, hdr_old, 'NAXIS'
  sxdelpar, hdr_old, 'NAXIS1'
  sxdelpar, hdr_old, 'NAXIS2'
  sxdelpar, hdr_old, 'BZERO'
  sxdelpar, hdr_old, 'BSCALE'
  sxdelpar, hdr_old, 'SIMPLE'
  sxdelpar, hdr_old, 'BITPIX'
  sxdelpar, hdr_old, 'EXTEND'
  hdr_out = [hdr1, hdr_old]
  sxaddpar, hdr_out, 'NEXP', nimgs
  sxaddpar, hdr_out, 'EXPTIME_TOT', exptime
  mwrfits, float(imgfinal), outfile2d, hdr_out, /create
  mwrfits, float(varfinal), outfile2d
  mwrfits, float(wavefinal), outfile2d
  mwrfits, objstruct_new, outfile2d

  RETURN
END
