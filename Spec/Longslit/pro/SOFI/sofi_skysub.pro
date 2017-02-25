FUNCTION SOFI_SKYSUB, afiles, bfiles, ivar = ivar $
                      , sky_resids = sky_resids, HDR = HDR $
                      , obj_pos = obj_pos, obj_neg = obj_neg $
                      , WAVEIMG = WAVEIMG $
                      , TELLURIC = TELLURIC1, CHK = CHK, WVCHK = WVCHK $
                      , VERBOSE = VERBOSE, targdir = targdir $
                      , AVG_SKY = AVG_SKY, SKY_MODEL = SKY_MODEL $
                      , NOMASK = NOMASK1 $
                      , SLITMASK = SLITMASK, tset_slits = tset_slits $
                      , wavefile = wavefile $
                      , peakthresh = peakthresh

;IF (N_params() LT 1) THEN BEGIN
;    print, 'Syntax: niri_skysub,filenames,skyfiles,flatfile,[OBJFIND=, CHK=,VERBOSE=,]'
;    return, 0
;endif

;flat = mrdfits(flatfile, 0)
;dims = size(flat, /dim)
slit_arcsec = 0.60d
plate_scale = 0.288d
slit = slit_arcsec/plate_scale        
FWHM = 3.5
sigma_psf = FWHM/2.35482D
;;pkwdth = slit
;;TOLER = slit/2.0D
nseq = n_elements(afiles)
fnseq = float(nseq)
fnseq2 = fnseq*fnseq

FOR ii = 0L, nseq-1L DO BEGIN
   sofi_proc, afiles[ii], aimg, ivar_a, hdr = hdr_a
   sofi_proc, bfiles[ii], bimg, ivar_b, hdr = hdr_b
   IF ii EQ 0 THEN BEGIN
      dims = size(aimg, /dim)
      nx = dims[0]
      ny = dims[1]
      a_stk = fltarr(nx*ny, nseq)
      ivar_a_stk = fltarr(nx*ny, nseq)
      mask_a_stk = fltarr(nx*ny, nseq)
      b_stk = fltarr(nx*ny, nseq)
      ivar_b_stk = fltarr(nx*ny, nseq)
      mask_b_stk = fltarr(nx*ny, nseq)
      hdr = hdr_a
   ENDIF
   a_stk[*, ii] = aimg
   ivar_a_stk[*, ii] = ivar_a
   mask_a_stk[*, ii] = (ivar_a LE 0.0)
   b_stk[*, ii] = bimg
   ivar_b_stk[*, ii] = ivar_b
   mask_b_stk[*, ii] = (ivar_b LE 0.0)
ENDFOR
;; Create an average sky for the wavelengths
skystack = fltarr(nx*ny, 2*nseq)
maskstack = fltarr(nx*ny, 2*nseq)
skystack[*, 0:(nseq-1L)] = a_stk[*, 0:nseq-1L]
skystack[*, nseq:(2*nseq-1L)] = b_stk[*, 0:nseq-1L]
maskstack[*, 0:(nseq-1L)] =  mask_a_stk[*, 0:nseq-1L]
maskstack[*, nseq:(2*nseq-1L)] = mask_b_stk[*, 0:nseq-1L]

avg_sky =  reform(djs_avsigclip(skystack, 2, sigrej = sigrej, inmask = maskstack, outmask = outmask), nx, ny)
smashmask = reform((total(outmask, 2) EQ 2*nseq), nx, ny) 
finalmask = (smashmask EQ 0)

;; ????? In the future make it so that the Telluric is read in here
;; and translated and used as slit boundaries. Then the telluric will
;; be used for tracing. 
tset_slits = niri_slitset(nx, ny)
slitmask = long_slits2mask(tset_slits)
ximg = long_slits2x(tset_slits)
y_img = replicate(1.0, nx)#findgen(ny)

CHK = 1
;;chk_trc = 1
IF KEYWORD_SET(TELLURIC1) THEN BEGIN
    TELLURIC = TELLURIC1
;;    NOMASK = 1
ENDIF ELSE TELLURIC = 0
IF NOT KEYWORD_SET(SKYBUFFER) THEN SKYBUFFER = 35L ;; 10" for SOFI
IF KEYWORD_SET(TELLURIC) THEN SKYBUFFER = 100L
;;chk_trc = 1

;; Create wavelength img 
qafile = targdir + '/waveQA-' + gnirs_fileprefix(afiles[0]) + '.ps'
IF KEYWORD_SET(TELLURIC) THEN BEGIN
   sofi_proc, telluric, wvimg, hdr = hdr_wv 
;              , pixflatfile = pixflatfile $
;              , illumflatfile = illumflatfile 
;; if we have a telluric std, use a science image specified by telluric
ENDIF ELSE BEGIN
   wvimg = avg_sky 
   hdr_wv = hdr_a
ENDELSE
waveimg = sofi_waveimg(wvimg, tset_slits, hdr_wv, piximg = piximg $
                       , qafile = qafile $
                       , CHK = CHK_TRC)
;;pixset = long_wavepix(avg_sky, tset_slits, FWHM = FWHM, pkwdth = pkwdth $
;;                      , toler = toler, chk = chk_trc)
;;piximg = long_wpix2image(pixset, tset_slits)
;IF NOT KEYWORD_SET(TELLURIC) THEN BEGIN
;   mwrfits, piximg, wavefile, /create
;   mwrfits, pixset, wavefile
;ENDIF
inslit = where(slitmask NE 0)
IF KEYWORD_SET(NOMASK1) THEN NOMASK = NOMASK1 $
ELSE NOMASK = 0


IF nseq GT 1 THEN aden = total((mask_a_stk EQ 0), 2) ELSE aden = (mask_a_stk EQ 0)
IF nseq GT 1 THEN a_bar = fnseq*(aden GT 0.0)*total(a_stk*(mask_a_stk EQ 0), 2)/(aden + (aden EQ 0.0)) $
ELSE a_bar = aden*a_stk
var_a_stk = (ivar_a_stk GT 0.0)/(ivar_a_stk + (ivar_a_stk LE 0.0))
IF nseq GT 1 THEN var_a_bar = $
   fnseq2*(aden GT 0.0)*total(var_a_stk*(mask_a_stk EQ 0), 2)/(aden + (aden EQ 0.0))^2 $
ELSE var_a_bar = aden*var_a_stk
ivar_a_bar = (var_a_bar GT 0.0)/(var_a_bar + (var_a_bar LE 0.0))
;;;
IF nseq GT 1 THEN bden = total((mask_b_stk EQ 0), 2) ELSE bden = (mask_b_stk EQ 0)
IF nseq GT 1 THEN b_bar = fnseq*(bden GT 0.0)*total(b_stk*(mask_b_stk EQ 0), 2)/(bden + (bden EQ 0.0)) $
ELSE b_bar = bden*b_stk
var_b_stk = (ivar_b_stk GT 0.0)/(ivar_b_stk + (ivar_b_stk LE 0.0))
IF nseq GT 1 THEN var_b_bar = $
   fnseq2*(bden GT 0.0)*total(var_b_stk*(mask_b_stk EQ 0), 2)/(bden + (bden EQ 0.0))^2 $
ELSE var_b_bar = bden*var_b_stk
ivar_b_bar = (var_b_bar GT 0.0)/(var_b_bar + (var_b_bar LE 0.0))

a_bar = reform(a_bar, nx, ny)
b_bar = reform(b_bar, nx, ny)

a_min_b = a_bar-b_bar
var = (ivar_a_bar GT 0.0)/(ivar_a_bar + (ivar_a_bar LE 0.0)) + $
      (ivar_b_bar GT 0.0)/(ivar_b_bar + (ivar_b_bar LE 0.0)) 
ivar = reform((var GT 0.0)/(var + (var LE 0.0)), nx, ny)

IF KEYWORD_SET(TELLURIC) THEN BEGIN
;   Skip second-pass for standard stars and bright objects.
    print, "    Bright calibration source; skipping 2nd pass..."
    sky_model = avg_sky
    sky_resids = 0.0*avg_sky
 ENDIF ELSE BEGIN
;      Otherwise do 2nd pass
    print, "  Second-pass: Global sky subtraction..."
    sky_resids = fltarr(nx, ny)
    buffer = 4
    inslit = where(slitmask NE 0 AND $
                   y_img GT 0 AND $ 
                   y_img LT 1024, nord)
    fitpix  = where(slitmask NE 0 AND $
                    finite(a_min_b) EQ 1 AND $
                    finite(ivar) EQ 1  AND $
                    ximg GT 0.05  AND $
                    ximg LT 0.95  AND $
                    abs(a_min_b) LE 5.4d4 AND $
                    ivar GT 0.0, npix)
    psort = sort(piximg[fitpix])
    sset = bspline_iterfit(piximg[fitpix[psort]], a_min_b[fitpix[psort]] $
                           , invvar = (ivar[fitpix[psort]] GT 0.0) $
                           , upper = 3, lower = 3 $
                           , bkspace = 1.1D, maxiter = 20, maxrej = 10)
    sky_resids[inslit] = bspline_valu(piximg[inslit], sset)
    IF KEYWORD_SET(CHK) THEN BEGIN
        plotx = piximg[fitpix[psort]]
        ploty = a_min_b[fitpix[psort]]
        rms = sqrt(djs_median(ploty^2))
        x_splot, plotx, ploty, psym1 = 3, ymnx = [-30.0*rms, 30.0*rms] $
                 , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
;, xrange = [min(plotx), max(plotx)] $

        wait, 1.5
    ENDIF
 ENDELSE

; now break the image up into two quadrants and subtract out any systematic
; horizontal features for each chip
;;a_min_b = a_min_
;;img_minsky = sciimg - avg_sky - global_sky
;smash_left = djs_avsigclip(img_minsky[0:(nx/2-1L), *], 1, sigrej = 2.0)
;corr_left = smash_left ## replicate(1.0, nx/2)
;smash_right = djs_avsigclip(img_minsky[nx/2:*, *], 1, sigrej = 2.0)
;corr_right = smash_right ## replicate(1.0, nx/2)
;corr = fltarr(nx, ny)
;corr[0:(nx/2-1), *] = corr_left
;corr[nx/2:*, *] = corr_right
;avg_sky = avg_sky + corr
;img_minsky = sciimg - avg_sky - global_sky
;img_minsky = img_minsky - corr

print, "Finding objects in sky-subtracted image"
; SOFI plate scale is 0.288. Assuming 1.0" seeing FWHM = 3.5 pix
a_min_b_fit = a_min_b - sky_resids
obj_pos = long_objfind(a_min_b_fit, tset_slits = tset_slits $
                       , FWHM = FWHM, OBJMASK = OBJMASK_POS $
                       , NPERSLIT = KEYWORD_SET(TELLURIC) $
                       , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                       , peakthresh = peakthresh)
IF NOT KEYWORD_SET(OBJ_POS) THEN BEGIN
   sky_model = avg_sky 
   return,  a_min_b 
ENDIF
fwhm = total(obj_pos.FWHM)/n_elements(obj_pos.FWHM)
; Search for negative objects that might be left over in the avg sky
obj_neg = long_objfind(-a_min_b_fit, tset_slits = tset_slits $
                       , FWHM = fwhm, OBJMASK = OBJMASK_NEG $
                       , NPERSLIT = KEYWORD_SET(TELLURIC) $
                       , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                       , peakthresh = peakthresh)
IF KEYWORD_SET(obj_neg) THEN objtot = [obj_pos, obj_neg] $
ELSE objtot = obj_pos
OBJMASK = (OBJMASK_POS EQ 1) OR (OBJMASK_NEG EQ 1)

;nobj_tot = n_elements(objtot)
nobj = n_elements(objtot)
IF NOT KEYWORD_SET(NOMASK) THEN BEGIN
   ;;local_sky = global_sky
    print, "Doing third pass local sky subtraction with object masking"
    FOR iobj = 0L, nobj-1L DO BEGIN
        skymask = lonarr(nx, ny)
        left  = objtot[iobj].xpos - SKYBUFFER 
        right = objtot[iobj].xpos + SKYBUFFER 
        FOR j = 0L, ny-1L DO BEGIN 
           xmin = floor(left[j]) >  0
           xmax = ceil(right[j]) < nx
           skymask[xmin:xmax, j] = 1L
           skymask[xmin:xmax, j] = 1L
        ENDFOR
        inobj = where(slitmask NE 0 AND   $
                      skymask  EQ 1 AND   $
                      y_img GT 0 AND      $ 
                      y_img LT 1024, nord)
        fitpix = where(slitmask NE 0 AND  $
                       skymask  EQ 1 AND  $
                       objmask  EQ 0 AND  $
                       finite(a_min_b) EQ 1 AND $
                       finite(ivar) EQ 1  AND $
                       ximg GT 0.05 AND $
                       ximg LT 0.95 AND $
                       abs(a_min_b) LE 5.4d4 AND $
                       ivar GT 0.0, npix)
        psort = sort(piximg[fitpix])
        ;; no weights on the sky-subtraction for now. 
        sset = bspline_iterfit(piximg[fitpix[psort]] $
                               , a_min_b[fitpix[psort]] $
                               , invvar = (ivar[fitpix[psort]] GT 0.0) $
                               , upper = 3, lower = 3 $
                               , bkspace = 1.1D, maxiter = 20, maxrej = 10)
        sky_resids[inobj] = bspline_valu(piximg[inobj], sset)
        IF KEYWORD_SET(CHK) THEN BEGIN
           plotx = piximg[fitpix[psort]]
           ploty = a_min_b[fitpix[psort]]
           rms = sqrt(djs_median(ploty^2))
           x_splot, plotx, ploty, psym1 = 3 $
                    , ymnx = [-10.0*rms, 10.0*rms] $
                    , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
;                     , xrange = [min(plotx), max(plotx)] $

        ENDIF
    ENDFOR
    sky_model = avg_sky
ENDIF ELSE BEGIN
   sky_model = avg_sky 
ENDELSE

print, "Finding objects in sky-subtracted image, second pass"
; SOFI plate scale is 0.288. Assuming 1.0" seeing FWHM = 3.5 pix
a_min_b_fit = a_min_b - sky_resids
obj_pos = long_objfind(a_min_b_fit, tset_slits = tset_slits $
                       , FWHM = FWHM, OBJMASK = OBJMASK_POS $
                       , NPERSLIT = KEYWORD_SET(TELLURIC) $
                       , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                       , peakthresh = peakthresh)
IF NOT KEYWORD_SET(OBJ_POS) THEN BEGIN
   sky_model = avg_sky 
   return,  a_min_b
ENDIF
fwhm = total(obj_pos.FWHM)/n_elements(obj_pos.FWHM)
; Search for negative objects that might be left over in the avg sky
obj_neg = long_objfind(-a_min_b_fit, tset_slits = tset_slits $
                       , FWHM = fwhm, OBJMASK = OBJMASK_NEG $
                       , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                       , NPERSLIT = KEYWORD_SET(TELLURIC) $
                       , peakthresh = peakthresh)
IF KEYWORD_SET(obj_neg) THEN objtot = [obj_pos, obj_neg] $
ELSE objtot = obj_pos
OBJMASK = (OBJMASK_POS EQ 1) OR (OBJMASK_NEG EQ 1)

;;waveimg = piximg
IF keyword_set(CHK) THEN BEGIN
   xatv, (a_min_b-sky_resids)*sqrt(ivar), wv = waveimg, /block
 ENDIF

RETURN, a_min_b
END
