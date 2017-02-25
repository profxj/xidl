FUNCTION NIRSPEC_SKYSUB, filename, skyfile,tset_slits $
                         , PIXFLATFILE=PIXFLATFILE $
                         , ILLUMFLATFILE=ILLUMFLATFILE,NOMASK = NOMASK1 $
                         , sciimg = sciimg, ivar_diff = ivar_diff $
                         , objstruct = obj_pos, WAVEIMG = WAVEIMG $
                         , SLITMASK = SLITMASK, MXSHFT=mxshft $
                         , TELLURIC = TELLURIC1, CHK = CHK, WVCHK = WVCHK $
                         , VERBOSE = VERBOSE, HDR = HDR, targdir = targdir $
                         , SKYIMG = SKYIMG, SIMPLE_SUB=simple_sub $
                         , FILTER=filter
  
IF NOT KEYWORD_SET(SKYBUFFER) THEN SKYBUFFER = 50L

IF KEYWORD_SET(NOMASK1) THEN NOMASK = NOMASK1 $
ELSE NOMASK = 0

IF KEYWORD_SET(TELLURIC1) THEN BEGIN
    TELLURIC = TELLURIC1
    NOMASK = 1
ENDIF ELSE TELLURIC = 0



;; The linear limit for NIRSPEC is 18,000 ADU. The gain is 5.0 so in
;; e- the linear limit is 90,000 e-. 
saturate = 3.0d5 ;; mask saturation in counts

nirspec_proc, filename, sciimg, sciivar, hdr = hdr, pixflatfile = pixflatfile $
              ,illumflatfile=illumflatfile
;; read in sky frame
nirspec_proc, skyfile, skyimg, skyivar, hdr = hdr1, pixflatfile = pixflatfile $
              ,illumflatfile=illumflatfile

dims = size(sciimg, /dimen)
nx = dims[0]
ny = dims[1]
;; Account for flexure between trace flat and data
IF NOT KEYWORD_SET(NOSHIFT) AND NOT KEYWORD_SET(TELLURIC) THEN $
   xshift = long_xcorr_slits(sciimg, tset_slits, /shift)
;; Generate left and right edge of slits
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge
slitmask = long_slits2mask(tset_slits)
ximg = long_slits2x(tset_slits)

AMP_SHIFT=10L
slitmask_left  = long_slits2mask(tset_slits, xshift = -amp_shift)
slitmask_right = long_slits2mask(tset_slits, xshift =  amp_shift)
ampmask = (slitmask_left EQ 0 AND slitmask_right EQ 0)
;; Correct for zero point amplifier offsets between quadrants
xval = findgen(nx) # replicate(1.0, ny)
yval = replicate(1.0, nx) # findgen(ny)
;;  quadrants
ul = xval GE 0    AND xval LE nx/2 AND yval GE ny/2 AND yval LE ny-1L
ur = xval GE nx/2 AND xval LE nx-1 AND yval GE ny/2 AND yval LE ny-1L
ll = xval GE 0    AND xval LE nx/2 AND yval GE 0    AND yval LE ny/2
lr = xval GE nx/2 AND xval LE nx-1 AND yval GE 0    AND yval LE ny/2
ul_ind = WHERE(ul AND ampmask)
ur_ind = WHERE(ur AND ampmask)
ll_ind = WHERE(ll AND ampmask)
lr_ind = WHERE(lr AND ampmask)
;; remove offset from sciimg
djs_iterstat, sciimg[ul_ind], mean = mean_ul, median = med_ul, sigrej = 2.0
sciimg[WHERE(ul)] -=  med_ul
djs_iterstat, sciimg[ur_ind], mean = mean_ur, median = med_ur, sigrej = 2.0
sciimg[WHERE(ur)] -=  med_ur 
djs_iterstat, sciimg[ll_ind], mean = mean_ll, median = med_ll, sigrej = 2.0
sciimg[WHERE(ll)] -=  med_ll 
djs_iterstat, sciimg[lr_ind], mean = mean_lr, median = med_lr, sigrej = 2.0
sciimg[WHERE(lr)] -=  med_lr
;; remove offset from skyimg
djs_iterstat, skyimg[ul_ind], mean = mean_ul, median = med_ul, sigrej = 2.0
skyimg[WHERE(ul)] -=  med_ul
djs_iterstat, skyimg[ur_ind], mean = mean_ur, median = med_ur, sigrej = 2.0
skyimg[WHERE(ur)] -=  med_ur 
djs_iterstat, skyimg[ll_ind], mean = mean_ll, median = med_ll, sigrej = 2.0
skyimg[WHERE(ll)] -=  med_ll 
djs_iterstat, skyimg[lr_ind], mean = mean_lr, median = med_lr, sigrej = 2.0
skyimg[WHERE(lr)] -=  med_lr

global_sky=fltarr(nx,ny)

;; Determine slit width from header
slit_hdr   =  strtrim(sxpar(hdr, 'SLITNAME'))
slit_str = strmid(slit_hdr,3,7)
slit_arcsec = double(slit_str[0])
slit = slit_arcsec/0.143D 
FWHM = slit
sigma_psf = FWHM/2.35482D
psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
FWHM_OBJ=1.2*slit ;; assume psf is 1.2 times slit width 

;; Create wavelength img 
qafile = targdir + '/waveQA-' + gnirs_fileprefix(filename[0]) + '.ps'
IF KEYWORD_SET(TELLURIC) THEN BEGIN
   nirspec_proc, telluric, wvimg, hdr = hdr_wv $
                 , pixflatfile = pixflatfile $
                 , illumflatfile=illumflatfile 
;; if we have a telluric std, use a science image specified by telluric
ENDIF ELSE BEGIN
   wvimg = sciimg 
   hdr_wv=hdr
ENDELSE

waveimg=nirspec_waveimg(wvimg,tset_slits,hdr_wv,piximg=piximg,qafile=qafile,MXSHFT=mxshft, $
                       FILTER=filter)

img_minsky = sciimg - skyimg
var_sci    = (sciivar GT 0.0)/(sciivar + (sciivar LE 0.0))
var_sky    = (skyivar GT 0.0)/(skyivar + (skyivar LE 0.0))
var_diff   = var_sci + var_sky
ivar_diff  = (var_diff GT 0.0)/(var_diff + (var_diff LE 0.0))
mask  = sciivar GT 0 AND skyivar GT 0 AND $
        sciimg LT saturate AND skyimg LT saturate 
ivar_diff=ivar_diff*mask

;; First pass, global sky subtraction
IF KEYWORD_SET(TELLURIC) THEN BEGIN
   ;;   Skip second-pass for standard stars and bright objects.
   print, "    Bright calibration source; skip B-spline sky subtraction..."
   sky_model = skyimg
   cr_sci = psf_reject_cr(img_minsky, ivar_diff, psfvals $
                          , nsigma=5.0, niter=10 $
                          , satmask = (img_minsky LT -10))
   cr_sky = psf_reject_cr(-img_minsky, ivar_diff, psfvals $
                          , nsigma=5.0, niter=10 $
                          , satmask = (img_minsky LT -10))
   icr = WHERE(cr_sci OR cr_sky,ncr)
   IF ncr GT 0 THEN mask[icr]=0
   ivar_diff=ivar_diff*mask
   obj_pos = long_objfind(mask*img_minsky, tset_slits = tset_slits $
                          , FWHM = FWHM_OBJ, OBJMASK = OBJMASK_POS $
                          , peakthresh = 0.1, /SILENT $
                          , OBJTHRESH = 0.1D)
ENDIF ELSE BEGIN
   print, "Doing global B-spline sky subtraction"
   inslit = where(slitmask NE 0)
   fitpix  = where(slitmask NE 0 AND $
                   finite(img_minsky) EQ 1 AND $
                   finite(sciivar) EQ 1  AND $
                   finite(skyivar) EQ 1  AND $
                   ximg GT 0.05  AND $
                   ximg LT 0.95  AND $
                   mask GT 0,npix)
   psort = sort(piximg[fitpix])
   sset = bspline_iterfit(piximg[fitpix[psort]], img_minsky[fitpix[psort]] $
                          , invvar = (mask[fitpix[psort]] GT 0) $
                          , upper = 3, lower = 3 $
                          , bkspace = 1.1D, maxiter = 20, maxrej = 10)
   global_sky[inslit]= bspline_valu(piximg[inslit], sset)
   img_minsky = sciimg - skyimg - global_sky
;; zap cosmics in both images
   cr_sci = psf_reject_cr(img_minsky, ivar_diff, psfvals $
                          , nsigma=5.0, niter=10 $
                          , satmask = (img_minsky LT -10))
   cr_sky = psf_reject_cr(-img_minsky, ivar_diff, psfvals $
                          , nsigma=5.0, niter=10 $
                          , satmask = (img_minsky LT -10))
   icr = WHERE(cr_sci OR cr_sky,ncr)
   IF ncr GT 0 THEN mask[icr]=0
   ivar_diff=ivar_diff*mask
;; Find and mask negative objects in the differenced frame
   obj_pos = long_objfind(mask*img_minsky, tset_slits = tset_slits $
                          , FWHM = FWHM_OBJ, OBJMASK = OBJMASK_POS $
                          , peakthresh = 0.1, /SILENT $
                          , OBJTHRESH = 0.1D, SIMPLE_SUB=simple_sub)
   obj_neg = long_objfind(-mask*img_minsky, tset_slits = tset_slits $
                          , FWHM = FWHM_OBJ, OBJMASK = OBJMASK_NEG $
                          , peakthresh = 0.1, /SILENT $
                          , OBJTHRESH = 0.1D, SIMPLE_SUB=simple_sub)
   OBJMASK = (OBJMASK_POS EQ 1) OR (OBJMASK_NEG EQ 1)
   ;stop
;; Second pass, global sky subtraction
;; Redo with cosmics and object traces masked
   img_minsky = sciimg - skyimg
   fitpix  = where(slitmask NE 0 AND $
                   finite(img_minsky) EQ 1 AND $
                   finite(sciivar) EQ 1  AND $
                   finite(skyivar) EQ 1  AND $
                   ximg GT 0.05  AND $
                   ximg LT 0.95  AND $
                   objmask EQ 0 AND $
                   mask GT 0, npix)
   psort = sort(piximg[fitpix])
   sset = bspline_iterfit(piximg[fitpix[psort]], img_minsky[fitpix[psort]] $
                          , invvar = (mask[fitpix[psort]] GT 0) $
                          , upper = 3, lower = 3 $
                          , bkspace = 1.1D, maxiter = 20, maxrej = 10)
   global_sky[inslit] = bspline_valu(piximg[inslit], sset)
   IF KEYWORD_SET(CHK) THEN BEGIN
      plotx = piximg[fitpix[psort]]
      ploty = img_minsky[fitpix[psort]]
      rms = sqrt(djs_median(ploty^2))
      x_splot, plotx, ploty, psym1 = 3, ymnx = [-30.0*rms, 30.0*rms] $
               , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
      wait, 1.5
   ENDIF
   img_minsky = sciimg - skyimg - global_sky
;; re-trace the positive object
   obj_pos = long_objfind(img_minsky, tset_slits = tset_slits $
                          , FWHM = FWHM_OBJ, OBJMASK = OBJMASK_POS $
                          , peakthresh = 0.1, /SILENT $
                          , OBJTHRESH = 0.1D, simple_sub=simple_sub)
   ;stop
ENDELSE
;; Third pass, local sky-subtraction
IF NOT KEYWORD_SET(NOMASK) THEN BEGIN
   img_minsky = sciimg - skyimg
   IF KEYWORD_SET(OBJ_POS) THEN nobj = n_elements(obj_pos) $
   ELSE nobj=0
   local_sky = global_sky
   print, "Doing B-spline local sky subtraction with object masking"
   FOR iobj = 0L, nobj-1L DO BEGIN
      skymask = lonarr(nx, ny)
      left  = obj_pos[iobj].xpos - SKYBUFFER 
      right = obj_pos[iobj].xpos + SKYBUFFER 
      FOR j = 0L, ny-1L DO BEGIN 
         xmin = floor(left[j]) >  0
         xmax = ceil(right[j]) < nx
         skymask[xmin:xmax, j] = 1L
         skymask[xmin:xmax, j] = 1L
      ENDFOR
      inobj = WHERE(slitmask GT 0 AND skymask EQ 1,nin)
      fitpix  = where(slitmask NE 0 AND $
                      finite(img_minsky) EQ 1 AND $
                      finite(sciivar) EQ 1  AND $
                      finite(skyivar) EQ 1  AND $
                      ximg GT 0.05  AND $
                      ximg LT 0.95  AND $
                      objmask EQ 0 AND $
                      skymask EQ 1 AND $
                      mask GT 0, npix)
      psort = sort(piximg[fitpix])
      ;; no weights on the sky-subtraction for now. 
      sset = bspline_iterfit(piximg[fitpix[psort]] $
                             , img_minsky[fitpix[psort]] $
                             , invvar = (mask[fitpix[psort]] GT 0) $
                             , upper = 3, lower = 3 $
                             , bkspace = 1.1D, maxiter = 20, maxrej = 10)
      local_sky[inobj] = bspline_valu(piximg[inobj], sset)
      IF KEYWORD_SET(CHK) THEN BEGIN
         plotx = piximg[fitpix[psort]]
         ploty = img_minsky[fitpix[psort]]
         rms = sqrt(djs_median(ploty^2))
         x_splot, plotx, ploty, psym1 = 3 $
                  , ymnx = [-10.0*rms, 10.0*rms] $
                  , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
      ENDIF
   ENDFOR
   sky_model = skyimg + local_sky
ENDIF 

IF keyword_set(CHK) THEN BEGIN
   xatv, (sciimg-sky_model)*sqrt(ivar_diff)*(slitmask GT 0) $
         , wv = waveimg, /block
ENDIF

RETURN, sky_model
END
