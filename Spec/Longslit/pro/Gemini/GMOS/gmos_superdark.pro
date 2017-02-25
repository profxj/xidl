PRO GMOS_SUPERDARK, biasfiles, darkfiles, outbias, outdark, sigrej = sigrej $
                    , IMAGING = IMAGING

  IF NOT KEYWORD_SET(SIGREJ) THEN SIGREJ = 3.0d
  ;; Make the superbias
  long_superbias, biasfiles, outbias
  ;; Make the superdark
  long_superbias, darkfiles, outdark
  
  superbias = mrdfits(outbias, 0)
  superdark = mrdfits(outdark, 0)
  
  dim = size(superbias, /dim)
  nx = dim[0]
  ny = dim[1]
  nxby3 = nx/3L
  xarr = lindgen(nx) # replicate(1.0, ny)
  yarr = replicate(1.0, nx) # lindgen(ny)
  
  IF      (ny-2L*9L)  MOD  512 EQ 0 THEN specbin = 4 $
  ELSE IF (ny-2L*18L) MOD  512 EQ 0 THEN specbin = 2 $
  ELSE IF (ny-2L*36L) MOD  512 EQ 0 THEN specbin = 1 $
  ELSE message, 'Transformation not supported for your binning'
  
  IF specbin EQ 4 THEN BEGIN
     ygap1 = 9L
     ygap2 = 9L
  ENDIF ELSE IF specbin EQ 2 THEN BEGIN
     ygap1 = 18L
     ygap2 = 18L
  ENDIF ELSE IF specbin EQ 1 THEN BEGIN
     ygap1 = 36L
     ygap2 = 36L
  ENDIF ELSE message, 'Your binning is not supported'
  
  IF NOT KEYWORD_SET(IMAGING) THEN nxmax = 2*nxby3-1L $
  ELSE nxmax = nx-1L
  ny_raw = (ny-ygap1-ygap2)/3
  chipmask = lonarr(nx, ny)
  chipmask[0:nxmax, 0:ny_raw-1L] = 3
  chipmask[0:nxmax, ny_raw + ygap1:2*ny_raw + ygap1 - 1L] = 2
  chipmask[0:nxmax, 2*ny_raw + ygap1 + ygap2:*] = 1
;;   Take the PSF width to be that of the spectral direction.  
;;   This prevents the routine from rejecting sky lines
;;   This is a description of the 3x3 core of the 2D PSF for reject_cr.pro
;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1] 
;;    PSFVALS[0]          1.   PSFVALS[0]
;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
  sigma_psf = 3.0d/2.35482D
  psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
  ncr = 0
  ivar  =  fltarr(nx, ny)
  mean  =  fltarr(nx, ny) 
  diff = superdark - superbias
;; Chip 1
  ichip1 = WHERE(chipmask EQ 1)
  djs_iterstat, diff[ichip1], med = med1, sigma = sigma1 
  mean[ichip1] = med1
  ivar[ichip1] = 1.0/sigma1^2
;; Chip 2
  ichip2 = WHERE(chipmask EQ 2)
  djs_iterstat, diff[ichip2], med = med2, sigma = sigma2
  mean[ichip2] = med2
  ivar[ichip2] = 1.0/sigma2^2
;; Chip 3
  ichip3 = WHERE(chipmask EQ 3)
  djs_iterstat, diff[ichip3], med = med3, sigma = sigma3
  mean[ichip3] = med3
  ivar[ichip3] = 1.0/sigma3^2
  
  med_width = 9.0
  darkmask = lonarr(nx, ny) + 1L
  niter = 3L
;; Mask convention is 1 = good, 0 = bad
  FOR iter = 0L, niter-1L DO BEGIN
     ;; sharp feature CR rejection
     cr_mask = psf_reject_cr(diff - mean, ivar*(chipmask EQ 0.0)*(ivar EQ 0.0) $
                             , psfvals, satmask = (diff-mean) GT 1d5 $
                             , nsig = sigrej)
     ;; threshold rejection
     sig = (ivar GT 0.0)/sqrt(ivar + (ivar EQ 0.0))
     sig_mask = (diff - mean) LT sigrej*sig
     ;; Define the mask
     darkmask = (cr_mask EQ 0) AND sig_mask
     sig_res = med_width/4.0
     nhalf =  long(sig_res)*4L
     xkern = findgen(2*nhalf+1)-nhalf
     kernel = gauss1(xkern, [0.0, sig_res, 1.0])
     
     mean_vec = djs_avsigclip(diff[0:nxmax, *], 1 $
                              , inmask = (darkmask[0:nxmax, *] EQ 0))
     mean_vec[0L:ny_raw-1L] = convol(djs_median(mean_vec[0L:ny_raw-1L] $
                                                , width = med_width $
                                                , boundary = 'reflect') $
                                     , kernel, /edge_truncate)
     mean_vec[ny_raw + ygap1:2L*ny_raw + ygap1-1L] = $
        convol(djs_median(mean_vec[ny_raw + ygap1:2L*ny_raw + ygap1-1L] $
                          , width = med_width $
                          , boundary = 'reflect'), kernel, /edge_truncate)
     mean_vec[2L*ny_raw + ygap1 + ygap2:*] = $
        convol(djs_median(mean_vec[2L*ny_raw + ygap1 + ygap2:*] $
                          , width = med_width $
                          , boundary = 'reflect'), kernel, /edge_truncate)
     mean = replicate(1.0, nx) # mean_vec
     ;; tweak the noise with this darkmask
     IF iter EQ 0 THEN BEGIN
        ;; Chip 1
        ichip1 = WHERE(chipmask EQ 1 AND darkmask)
        djs_iterstat, diff[ichip1] - mean[ichip1], med = med1, sigma = sigma1 
        ivar[ichip1] = 1.0/sigma1^2
        ;; Chip 2
        ichip2 = WHERE(chipmask EQ 2 AND darkmask)
        djs_iterstat, diff[ichip2]-mean[ichip2], med = med2, sigma = sigma2
        ivar[ichip2] = 1.0/sigma2^2
        ;; Chip 3
        ichip3 = WHERE(chipmask EQ 3 AND darkmask)
        djs_iterstat, diff[ichip3]-mean[ichip3], med = med3, sigma = sigma3
        ivar[ichip3] = 1.0/sigma3^2
     ENDIF
  ENDFOR
  ;; set the mask to 1 in the 
  IF NOT KEYWORD_SET(IMAGING) THEN darkmask[2*nxby3:*, *] = 1L

  dark = mrdfits(outdark, 0, hdr)
  mwrfits, dark, outdark, hdr, /create
  mwrfits, darkmask, outdark
  
  RETURN
END
