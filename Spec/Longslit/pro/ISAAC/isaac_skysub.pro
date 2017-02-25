FUNCTION ISAAC_SKYSUB, filename, skyfile, scifile, darkfile = darkfile $
                       , tset_slits = tset_slit $
                       , PIXFLATFILE = PIXFLATFILE $
                       , ILLUMFLATFILE = ILLUMFLATFILE, NOMASK = NOMASK1 $
                       , sciimg = sciimg, ivar_diff = ivar_diff $
                       , objstruct = obj_pos, WAVEIMG = WAVEIMG $
                       , SLITMASK = SLITMASK $
                       , TELLURIC = TELLURIC1, USE_TELL_WAVE = USE_TELL_WAVE $
                       , CHK = CHK, WVCHK = WVCHK $
                       , VERBOSE = VERBOSE, HDR = HDR, targdir = targdir $
                       , SKYIMG = SKYIMG, LOCAL_SKY = LOCAL_SKY $
                       , STDTRACE = STDTRACE, CALIB = CALIB

IF NOT KEYWORD_SET(SKYBUFFER) THEN SKYBUFFER = 68L ;; 10" for ISAAC

IF KEYWORD_SET(NOMASK1) THEN NOMASK = NOMASK1 $
ELSE NOMASK = 0

IF KEYWORD_SET(TELLURIC1) THEN BEGIN
    TELLURIC = TELLURIC1
    NOMASK = 1
    SLITFILE = TELLURIC
ENDIF ELSE BEGIN 
   TELLURIC = 0
   SLITFILE = FILENAME
ENDELSE
bkspace = 1.1


isaac_slitmask, slitfile, tset_slits = tset_slits, slitmask = slitmask $
                , darkfile = darkfile

saturate = 4d4*4.5d ;; manual says 40,000 ADU and gain is 4.5 e^-/ADU 

isaac_proc, filename, sciimg, sciivar, hdr = hdr, pixflatfile = pixflatfile $
            , illumflatfile = illumflatfile, darkfile = darkfile
;; read in sky frame
isaac_proc, skyfile, skyimg, skyivar, hdr = hdr1, pixflatfile = pixflatfile $
            , illumflatfile = illumflatfile, darkfile = darkfile

dims = size(sciimg, /dimen)
nx = dims[0]
ny = dims[1]

NOSHIFT = 1
;; Account for flexure between trace flat and data
;IF NOT KEYWORD_SET(NOSHIFT) AND NOT KEYWORD_SET(TELLURIC) THEN $
;   xshift = long_xcorr_slits(sciimg, tset_slits, /shift)
;; Generate left and right edge of slits
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge
;;slitmask = long_slits2mask(tset_slits)
ximg = long_slits2x(tset_slits)

dim_slit = size(tset_slits[0].coeff, /dim)
IF n_elements(dim_slit) EQ 1 THEN nslit = 1 ELSE nslit = dim_slit[1]

global_sky = fltarr(nx, ny)

;; Determine slit width from header
slit_str = strcompress(esopar(hdr, 'HIERARCH ESO INS OPTI1 ID'), /rem)
CASE slit_str OF 
   'slit_1': slit_width = 1.0d
   'slit_0.6_tilted': slit_width = 0.6d
   ELSE: message, 'Unrecognized slit'
ENDCASE
plate_scale = 0.147d
fnslit = slit_width/plate_scale
FWHM = fnslit
sigma_psf = FWHM/2.35482D
psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
FWHM_OBJ = 1.2*fnslit ;; assume psf is 1.2 times slit width 

;; Create wavelength img 
;;qafile = targdir + '/waveQA-' + gnirs_fileprefix(filename[0]) + '.ps'
IF KEYWORD_SET(TELLURIC) THEN BEGIN 
   ;; if TELLURIC = 'use_tell_wave' then use the telluric for
   ;; wavelength calibration
   IF TELLURIC[0] EQ 'use_tell_wave' THEN BEGIN 
      nfiles = 2 
      use_tell_wave = 1
   ENDIF ELSE nfiles = n_elements(TELLURIC)
ENDIF ELSE nfiles = 2
skystack  = fltarr(nx*ny, nfiles)
var_stack = fltarr(nx*ny, nfiles)
maskstack = lonarr(nx*ny, nfiles)
IF KEYWORD_SET(TELLURIC) AND NOT KEYWORD_SET(USE_TELL_WAVE) THEN BEGIN
   ;; Read in sky files (this is wavelength files for TELLURIC)
   FOR j = 0L, nfiles-1L DO BEGIN
      isaac_proc, TELLURIC[j], imag1, ivar1, hdr = hdr1 $
                  , pixflatfile = pixflatfile $
                  , illumflatfile = illumflatfile, darkfile = darkfile
      skystack[*, j] = imag1
      var_stack[*, j] = float(ivar1 GT 0.0)/(ivar1 + (ivar1 LE 0.0))
      maskstack[*, j] = (ivar1 LE 0.0) ; opposite convention for avsigclip
   ENDFOR
   hdr_wv = hdr1
ENDIF ELSE BEGIN
   skystack[*, 0] = sciimg
   skystack[*, 1] = skyimg
   var_stack[*, 0] = float(sciivar GT 0.0)/(sciivar + (sciivar LE 0.0))
   var_stack[*, 1] = float(skyivar GT 0.0)/(skyivar + (skyivar LE 0.0))
   maskstack[*, 0] = (sciivar LE 0.0) ; opposite convention for avsigclip
   maskstack[*, 1] = (skyivar LE 0.0) ; opposite convention for avsigclip
   hdr_wv = hdr
ENDELSE
;; This sigrej is only for combining the wavelength solution images
IF (NOT keyword_set(sigrej1)) then begin
   if (nfiles LE 2) then sigrej = 1.0 $ 
;; Irrelevant for only 1 or 2 files
   else if (nfiles EQ 3) then sigrej = 1.1 $
   else if (nfiles EQ 4) then sigrej = 1.3 $
   else if (nfiles EQ 5) then sigrej = 1.6 $
   else if (nfiles EQ 6) then sigrej = 1.9 $
   else sigrej = 2.0
ENDIF ELSE sigrej = sigrej1
IF nfiles GT 1 THEN BEGIN
   wv_img = reform(djs_avsigclip(skystack, 2, sigrej = sigrej, inmask = maskstack $
                                 , outmask = outmask), nx, ny)
   nused = total(outmask EQ 0, 2)
   wv_var = total(var_stack*(outmask EQ 0), 2)/(nused^2 + (nused EQ 0))
   maskfinal = (total(outmask, 2) NE nfiles)
   wv_ivar = reform(float(wv_var GT 0.0)/(wv_var + (wv_var LE 0.0)), nx, ny)
ENDIF ELSE BEGIN
   wv_img = reform(skystack, nx, ny)
   wv_ivar = reform(float(var_stack GT 0.0)/(var_stack + (var_stack LE 0.0)) $
                    , nx, ny)
ENDELSE
   ;isaac_proc, telluric, wv_img, wv_ivar, hdr = hdr_wv $
   ;            , pixflatfile = pixflatfile $
   ;            , illumflatfile = illumflatfile 
;; if we have a telluric std, use a science image specified by the
;; telluric variable
waveimg = isaac_waveimg(wv_img, wv_ivar, tset_slits, hdr_wv, scifile $
                        , piximg = piximg, CHK = WVCHK, CALIB = CALIB)

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
                          , OBJTHRESH = 0.1D, STDTRACE = STDTRACE)
ENDIF ELSE BEGIN
   FOR islit = 0L, nslit-1L DO BEGIN
      print, "Doing global B-spline sky subtraction"
      inslit = where(slitmask EQ (islit +1L))
      fitpix  = where(slitmask EQ (islit + 1L) AND $
                      finite(img_minsky) EQ 1 AND $
                      finite(sciivar) EQ 1  AND $
                      finite(skyivar) EQ 1  AND $
                      ximg GT 0.05  AND $
                      ximg LT 0.95  AND $
                      mask GT 0, npix)
      psort = sort(piximg[fitpix])
      sset = bspline_iterfit(piximg[fitpix[psort]], img_minsky[fitpix[psort]] $
                             , invvar = (mask[fitpix[psort]] GT 0) $
                             , upper = 3, lower = 3 $
                             , bkspace = bkspace, maxiter = 20, maxrej = 10)
      global_sky[inslit] = bspline_valu(piximg[inslit], sset)
   ENDFOR
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
                          , OBJTHRESH = 0.1D, STDTRACE = STDTRACE)
   obj_neg = long_objfind(-mask*img_minsky, tset_slits = tset_slits $
                          , FWHM = FWHM_OBJ, OBJMASK = OBJMASK_NEG $
                          , peakthresh = 0.1, /SILENT $
                          , OBJTHRESH = 0.1D, STDTRACE = STDTRACE)
   OBJMASK = (OBJMASK_POS EQ 1) OR (OBJMASK_NEG EQ 1)
;; Second pass, global sky subtraction
;; Redo with cosmics and object traces masked
   img_minsky = sciimg - skyimg
   FOR islit = 0L, nslit-1L DO BEGIN
      inslit = where(slitmask EQ (islit +1L))
      fitpix  = where(slitmask EQ (islit + 1L) AND $
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
                             , bkspace = bkspace, maxiter = 20, maxrej = 10)
      global_sky[inslit] = bspline_valu(piximg[inslit], sset)
      IF KEYWORD_SET(CHK) THEN BEGIN
         plotx = piximg[fitpix[psort]]
         ploty = img_minsky[fitpix[psort]]
         rms = sqrt(djs_median(ploty^2))
         x_splot, plotx, ploty, psym1 = 3, ymnx = [-30.0*rms, 30.0*rms] $
                  , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
         wait, 1.5
      ENDIF
   ENDFOR
   img_minsky = sciimg - skyimg - global_sky
;; re-trace the positive object
   obj_pos = long_objfind(img_minsky, tset_slits = tset_slits $
                          , FWHM = FWHM_OBJ, OBJMASK = OBJMASK_POS $
                          , peakthresh = 0.1, /SILENT $
                          , OBJTHRESH = 0.1D, STDTRACE = STDTRACE)
   OBJMASK = (OBJMASK_POS EQ 1) OR (OBJMASK_NEG EQ 1)
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
         xmax = ceil(right[j]) < (nx-1L)
         skymask[xmin:xmax, j] = 1L
         skymask[xmin:xmax, j] = 1L
      ENDFOR
      inobj = WHERE(slitmask EQ obj_pos[iobj].SLITID AND skymask EQ 1,nin)
      fitpix  = where(slitmask EQ obj_pos[iobj].SLITID AND $
                      finite(img_minsky) EQ 1 AND $
                      finite(sciivar) EQ 1 AND $
                      finite(skyivar) EQ 1 AND $
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
                             , bkspace = bkspace, maxiter = 20, maxrej = 10)
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
