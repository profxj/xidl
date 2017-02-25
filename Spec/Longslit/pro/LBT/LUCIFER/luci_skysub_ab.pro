FUNCTION LUCI_SKYSUB_AB, a_min_b, ivar, tset_slits, piximg $
                         , STDTRACE = STDTRACE $
                         , obj_pos = obj_pos, obj_neg = obj_neg $
                         , FWHM = FWHM1 $
                         , OBJMASK_OLD = OBJMASK_OLD $
                         , TELLURIC = TELLURIC1, CHK = CHK $
                         , VERBOSE = VERBOSE $
                         , NOMASK = NOMASK1 $
                         , peakthresh = peakthresh, SOFI = SOFI,  SIZE_OBJMASK = SIZE_OBJMASK1

;IF (N_params() LT 1) THEN BEGIN
;    print, 'Syntax: niri_skysub,filenames,skyfiles,flatfile,[OBJFIND=, CHK=,VERBOSE=,]'
;    return, 0
;endif
  sigrej = 3
;flat = mrdfits(flatfile, 0)
;dims = size(flat, /dim)
IF KEYWORD_SET(SIZE_OBJMASK1) THEN SIZE_OBJMASK = SIZE_OBJMASK1 ELSE SIZE_OBJMASK = 1.5
bsp = 0.7D
slit_arcsec = 1.0d
IF KEYWORD_SET(SOFI) THEN plate_scale = 0.288d ELSE plate_scale = 0.250d
slit = slit_arcsec/plate_scale        
IF KEYWORD_SET(FWHM1) THEN FWHM = FWHM1 ELSE FWHM = 4.5
sigma_psf = FWHM/2.35482D
;;pkwdth = slit
;;TOLER = slit/2.0D
nseq = n_elements(afiles)
fnseq = float(nseq)
fnseq2 = fnseq*fnseq
diff_max = 2d5 ;;5.4d4

;; Maximum FWHM = 1.3" seeing
IF NOT KEYWORD_SET(MAXFWHM) THEN MAXFWHM = 1.3D/plate_scale
;; Minimum FWHM is 1 pixel
IF NOT KEYWORD_SET(MINFWHM) THEN MINFWHM = 1.0D


;; ------
;; Expand slit set to get left and right edge
if not keyword_set(ARCTRC_POS) then arctrc_pos = 0.5
dims = tset_slits[0].DIMS
nx = dims[0]
ny = dims[1]
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge
trace = (left_edge + right_edge) * arctrc_pos

nsamp = 50
maxflex = 5L
step = lindgen(2*MAXFLEX*nsamp) - MAXFLEX*nsamp 


;; ????? In the future make it so that the Telluric is read in here
;; and translated and used as slit boundaries. Then the telluric will
;; be used for tracing. 
slitmask = long_slits2mask(tset_slits)
ximg = long_slits2x(tset_slits)
y_img = replicate(1.0, nx)#findgen(ny)

;;CHK = 1
;;chk_trc = 1
IF KEYWORD_SET(TELLURIC1) THEN BEGIN
    TELLURIC = TELLURIC1
    NOMASK = 1
ENDIF ELSE TELLURIC = 0
IF NOT KEYWORD_SET(SKYBUFFER) THEN SKYBUFFER = 50L ;; 12" for LUCI
IF KEYWORD_SET(TELLURIC) THEN SKYBUFFER = 100L
;;chk_trc = 1


sky_resids = fltarr(nx, ny)
inslit = where(slitmask NE 0)
IF KEYWORD_SET(NOMASK1) THEN NOMASK = NOMASK1 $
ELSE NOMASK = 0

IF NOT KEYWORD_SET(TELLURIC) THEN BEGIN
;;  First pass sky-subtraction 
    print, "  First-pass: Global sky subtraction..."
    buffer = 4
    inslit = where(slitmask NE 0 AND $
                   y_img GT 0 AND $ 
                   y_img LT ny, nord)
    fitpix  = where(slitmask NE 0 AND $
                    finite(a_min_b) EQ 1 AND $
                    finite(ivar) EQ 1  AND $
                    ximg GT 0.05  AND $
                    ximg LT 0.95  AND $
                    abs(a_min_b) LE diff_max AND $
                    ivar GT 0.0, npix)
    psort = sort(piximg[fitpix])
    sset = bspline_iterfit(piximg[fitpix[psort]], a_min_b[fitpix[psort]] $
                           , invvar = (ivar[fitpix[psort]] GT 0.0) $
                           , upper = 3, lower = 3 $
                           , bkspace = bsp, maxiter = 20, maxrej = 10)
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
 ENDIF

print, "Finding objects in sky-subtracted image"
; SOFI plate scale is 0.288. Assuming 1.0" seeing FWHM = 3.5 pix
a_min_b_fit = a_min_b - sky_resids
IF NOT KEYWORD_SET(OBJ_POS) THEN $
   obj_pos = long_objfind(a_min_b_fit, tset_slits = tset_slits $
                          , FWHM = FWHM $
;;, OBJMASK = OBJMASK_POS $
                          , NPERSLIT = KEYWORD_SET(TELLURIC) $
                          , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                          , peakthresh = peakthresh, STDTRACE = STDTRACE)
IF KEYWORD_SET(OBJ_POS) THEN fwhm = total(obj_pos.FWHM)/n_elements(obj_pos.FWHM)
IF NOT KEYWORD_SET(OBJ_NEG) THEN $
   obj_neg = long_objfind(-a_min_b_fit, tset_slits = tset_slits $
                          , FWHM = fwhm $
;;, OBJMASK = OBJMASK_NEG $
                          , NPERSLIT = KEYWORD_SET(TELLURIC) $
                          , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                          , peakthresh = peakthresh, STDTRACE = STDTRACE)
IF KEYWORD_SET(obj_pos) AND KEYWORD_SET(OBJ_neg) THEN objtot = [obj_pos, obj_neg] $
ELSE IF KEYWORD_SET(OBJ_POS) AND NOT KEYWORD_SET(OBJ_NEG) THEN objtot = obj_pos $
ELSE IF NOT KEYWORD_SET(OBJ_POS) AND KEYWORD_SET(OBJ_NEG) THEN objtot = obj_neg $
ELSE return, sky_resids
;;IF NOT KEYWORD_SET(OBJMASK_OLD) THEN OBJMASK_OLD = (OBJMASK_POS EQ 1) OR (OBJMASK_NEG EQ 1)
 
nobj = n_elements(objtot)

IF NOT KEYWORD_SET(NOMASK) AND NOT KEYWORD_SET(TELLURIC) THEN BEGIN
    ;;local_sky = global_sky
    print, "Doing second pass local sky subtraction with object masking"
    ;; Creat my own object mask
    obj_fwhm = (objtot.FWHM > MINFWHM) < MAXFWHM
    left_obj  = objtot.xpos -  $
                SIZE_OBJMASK*obj_fwhm ## replicate(1.0D, ny)
    right_obj = objtot.xpos + $
                SIZE_OBJMASK*obj_fwhm ## replicate(1.0D, ny)
    objmask = lonarr(nx, ny)
    objmask[where(slitmask)] = 1L
    FOR iobj = 0L, nobj-1L DO BEGIN
       FOR j = 0, ny-1L DO BEGIN
          xmin_obj = floor(left_obj[j, iobj]) >  0
          xmax_obj = ceil(right_obj[j, iobj]) < (nx-1L)
          objmask[xmin_obj:xmax_obj, j] = 0L
       ENDFOR
    ENDFOR
    FOR iobj = 0L, nobj-1L DO BEGIN
       skymask = lonarr(nx, ny)
       left  = objtot[iobj].xpos - SKYBUFFER 
       right = objtot[iobj].xpos + SKYBUFFER 
       FOR j = 0L, ny-1L DO BEGIN 
          xmin = floor(left[j]) >  0
          xmax = ceil(right[j]) < (nx-1L)
          skymask[xmin:xmax, j] = 1L
          skymask[xmin:xmax, j] = 1L
       ENDFOR
       inobj = where(slitmask NE 0 AND   $
                     skymask  EQ 1 AND   $
                     y_img GT 0 AND      $ 
                     y_img LT ny, nord)
       fitpix = where(slitmask NE 0 AND  $
                      skymask  EQ 1 AND  $
                      objmask  EQ 1 AND  $
                      finite(a_min_b) EQ 1 AND $
                      finite(ivar) EQ 1  AND $
                      ximg GT 0.05 AND $
                      ximg LT 0.95 AND $
                      abs(a_min_b) LE diff_max AND $
                      ivar GT 0.0, npix)
       psort = sort(piximg[fitpix])
       ;; no weights on the sky-subtraction for now. 
       sset = bspline_iterfit(piximg[fitpix[psort]] $
                              , a_min_b[fitpix[psort]] $
                              , invvar = (ivar[fitpix[psort]] GT 0.0) $
                              , upper = 3, lower = 3 $
                              , bkspace = bsp, maxiter = 20, maxrej = 10)
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
 ENDIF
 
 IF keyword_set(CHK) THEN BEGIN
    xatv, (a_min_b-sky_resids)*sqrt(ivar), /block
 ENDIF

 RETURN, sky_resids
END

;print, "Finding objects in sky-subtracted image, second pass"
; SOFI plate scale is 0.288. Assuming 1.0" seeing FWHM = 3.5 pix
;a_min_b_fit = a_min_b - sky_resids
;IF NOT KEYWORD_SET(OBJ_POS) THEN $
;   obj_pos = long_objfind(a_min_b_fit, tset_slits = tset_slits $
;                         , FWHM = FWHM, OBJMASK = OBJMASK_POS $
;                          , NPERSLIT = KEYWORD_SET(TELLURIC) $
;                          , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
;                          , peakthresh = peakthresh)
;IF KEYWORD_SET(OBJ_POS) THEN fwhm = total(obj_pos.FWHM)/n_elements(obj_pos.FWHM)
; Search for negative objects that might be left over in the avg sky
;IF NOT KEYWORD_SET(OBJ_NEG) THEN $
;   obj_neg = long_objfind(-a_min_b_fit, tset_slits = tset_slits $
;                          , FWHM = fwhm, OBJMASK = OBJMASK_NEG $
;                          , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
;                          , NPERSLIT = KEYWORD_SET(TELLURIC) $
;                          , peakthresh = peakthresh)
;IF KEYWORD_SET(obj_neg) THEN objtot = [obj_pos, obj_neg] $
;ELSE objtot = obj_pos
;OBJMASK = (OBJMASK_POS EQ 1) OR (OBJMASK_NEG EQ 1)
