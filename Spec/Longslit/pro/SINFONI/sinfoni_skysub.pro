FUNCTION SINFONI_SKYSUB, objfile, skyfile, tset_slits, piximg $
                         , ivar = ivar $
                         , sky_resids = sky_resids $
                         , sky_model = sky_model $
                         , plate_scale = plate_scale $
                         , hdr_obj = hdr_obj $
                         , pixflatfile = pixflatfile $
                         , illumflatfile = illumflatfile $
                         , obj = obj_pos, obj_neg = obj_neg $
                         , FWHM = FWHM1 $
                         , OBJMASK = OBJMASK $
                         , TELLURIC = TELLURIC, CHK = CHK $
                         , VERBOSE = VERBOSE $
                         , NOMASK = NOMASK1 $
                         , peakthresh = peakthresh

  sigrej = 3
  bsp = 0.7D
;; Do first sky-subtraction, proc images but don't dark subtract
  sinfoni_proc, objfile, obj_img, ivar_obj, hdr = hdr_obj $
                , pixflat = pixflatfile, illumflat = illumflatfile
  sinfoni_proc, skyfile, sky_img, ivar_sky, hdr = hdr_sky $
                , pixflat = pixflatfile, illumflat = illumflatfile
  mask_obj = (ivar_obj GT 0.0)
  mask_sky = (ivar_sky GT 0.0)
  var_obj = (ivar_obj GT 0.0)/(ivar_obj + (ivar_obj LE 0.0))
  var_sky = (ivar_sky GT 0.0)/(ivar_sky + (ivar_sky LE 0.0))
  obj_min_sky = obj_img - sky_img
  var_obj_min_sky = var_obj + var_sky
  mask = (mask_obj GT 0.0) AND (mask_sky GT 0.0)
  ivar = float(mask)/(var_obj_min_sky + (var_obj_min_sky EQ 0.0))

  optic = strcompress(esopar(hdr_obj, 'HIERARCH ESO INS OPTI1 NAME'), /rem)
  grating = strcompress(esopar(hdr_obj, 'HIERARCH ESO INS GRAT1 NAME'), /rem)
  slit_width = double(optic)
  plate_scale = slit_width/2.0d

  dims = tset_slits[0].DIMS
  nx = dims[0]
  ny = dims[1]
   
  dimt = size(tset_slits.coeff, /dimen)
  nslit = dimt[1]
  
  
  sky_resids = fltarr(nx, ny)
  sky_model = sky_img
  ;; No residual subtraction for LGS mode or TELLURIC standards
  IF optic EQ '0.025' OR KEYWORD_SET(TELLURIC) THEN RETURN, obj_min_sky
  
  IF KEYWORD_SET(FWHM1) THEN FWHM = FWHM1 ELSE FWHM = 4.5
  diff_max = 2d5 ;;5.4d4
  IF NOT KEYWORD_SET(MAXFWHM) THEN MAXFWHM = 1.3D/plate_scale
  ;; Minimum FWHM is 1 pixel
  IF NOT KEYWORD_SET(MINFWHM) THEN MINFWHM = 1.0D
  

;; ------
;; Expand slit set to get left and right edge
  ;;if not keyword_set(ARCTRC_POS) then arctrc_pos = 0.5
   ;traceset2xy, tset_slits[0], rows, left_edge
   ;traceset2xy, tset_slits[1], rows, right_edge
   ;;trace = (left_edge + right_edge) * arctrc_pos
   
   ;;nsamp = 50
   ;;maxflex = 5L
   ;;step = lindgen(2*MAXFLEX*nsamp) - MAXFLEX*nsamp 


   slitmask = long_slits2mask(tset_slits)
   ximg = long_slits2x(tset_slits, edgmask = edgmask)
   y_img = replicate(1.0, nx)#findgen(ny)

;;CHK = 1
;;chk_trc = 1
;IF KEYWORD_SET(TELLURIC1) THEN BEGIN
;    TELLURIC = TELLURIC1
;    NOMASK = 1
;ENDIF ELSE TELLURIC = 0


IF KEYWORD_SET(NOMASK1) THEN NOMASK = NOMASK1 $
ELSE NOMASK = 0


IF NOT KEYWORD_SET(TELLURIC) THEN BEGIN
;;  First pass sky-subtraction
   print, "  First-pass: Global sky subtraction..."
   FOR islit = 0L, nslit-1L DO BEGIN 
      buffer = 4
      inslit = where(slitmask EQ (islit +1L), nin)
      fitpix  = where(slitmask EQ (islit + 1L) AND $
                      finite(obj_min_sky) EQ 1 AND $
                      finite(ivar) EQ 1  AND $
                      edgmask EQ 0 AND $
                      ;;ximg GT 0.05  AND $
                      ;;ximg LT 0.95  AND $
                      abs(obj_min_sky) LE diff_max AND $
                      ivar GT 0.0, npix)
      psort = sort(piximg[fitpix])
      sset = bspline_iterfit(piximg[fitpix[psort]], obj_min_sky[fitpix[psort]] $
                             , invvar = (ivar[fitpix[psort]] GT 0.0) $
                             , upper = 3, lower = 3 $
                             , bkspace = bsp, maxiter = 20, maxrej = 10)
      sky_resids[inslit] = bspline_valu(piximg[inslit], sset)
      IF KEYWORD_SET(CHK) THEN BEGIN
         plotx = piximg[fitpix[psort]]
         ploty = obj_min_sky[fitpix[psort]]
         rms = sqrt(djs_median(ploty^2))
         x_splot, plotx, ploty, psym1 = 3, ymnx = [-30.0*rms, 30.0*rms] $
                  , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
;, xrange = [min(plotx), max(plotx)] $
         ;;wait, 1.5
      ENDIF
   ENDFOR
ENDIF

print, "Finding objects in sky-subtracted image"
obj_min_sky_fit = obj_min_sky - sky_resids
IF NOT KEYWORD_SET(OBJ_POS) THEN $
   obj_pos = long_objfind(obj_min_sky_fit, tset_slits = tset_slits $
                          , FWHM = FWHM, OBJMASK = OBJMASK_POS $
                          , NPERSLIT = KEYWORD_SET(TELLURIC) $
                          , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                          , peakthresh = peakthresh)
IF KEYWORD_SET(OBJ_POS) THEN fwhm = total(obj_pos.FWHM)/n_elements(obj_pos.FWHM)
IF NOT KEYWORD_SET(OBJ_NEG) THEN $
   obj_neg = long_objfind(-obj_min_sky_fit, tset_slits = tset_slits $
                          , FWHM = fwhm, OBJMASK = OBJMASK_NEG $
                          , NPERSLIT = KEYWORD_SET(TELLURIC) $
                          , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                          , peakthresh = peakthresh)
IF KEYWORD_SET(obj_pos) AND KEYWORD_SET(OBJ_neg) THEN objtot = [obj_pos, obj_neg] $
ELSE IF KEYWORD_SET(OBJ_POS) AND NOT KEYWORD_SET(OBJ_NEG) THEN objtot = obj_pos $
ELSE IF NOT KEYWORD_SET(OBJ_POS) AND KEYWORD_SET(OBJ_NEG) THEN objtot = obj_neg $
ELSE return, sky_resids
IF NOT KEYWORD_SET(OBJMASK) THEN OBJMASK = (OBJMASK_POS EQ 1) OR (OBJMASK_NEG EQ 1)

 nobj = n_elements(objtot)
 IF NOT KEYWORD_SET(NOMASK) THEN BEGIN
    print, "Doing second pass local sky subtraction with object masking"
    ;; Construct the mask
    obj_fwhm = (objtot.FWHM > MINFWHM) < MAXFWHM
    left_obj  = objtot.xpos -  $
                1.3*obj_fwhm ## replicate(1.0D, ny)
    right_obj = objtot.xpos + $
                1.3*obj_fwhm ## replicate(1.0D, ny)
    skymask = lonarr(nx, ny)
    skymask[WHERE(slitmask)] = 1L
    ;; Loop over the objects to create skymask
    FOR iobj = 0L, nobj-1L DO BEGIN
       FOR j = 0, ny-1L DO BEGIN
          xmin_obj = floor(left_obj[j, iobj]) >  0
          xmax_obj = ceil(right_obj[j, iobj]) < (nx-1L)
          skymask[xmin_obj:xmax_obj, j] = 0L
       ENDFOR
    ENDFOR
    FOR islit = 0L, nslit-1L DO BEGIN 
       inslit = where(slitmask EQ (islit +1L), nin) ;;AND $
       fitpix  = where(slitmask EQ (islit + 1L) AND $
                       finite(obj_min_sky) EQ 1 AND $
                       finite(ivar) EQ 1  AND $
                       edgmask EQ 0 AND $
                       abs(obj_min_sky) LE diff_max AND $
                       ivar GT 0.0 AND $
                       skymask EQ 1, npix)
       psort = sort(piximg[fitpix])
       sset = bspline_iterfit(piximg[fitpix[psort]], obj_min_sky[fitpix[psort]] $
                              , invvar = (ivar[fitpix[psort]] GT 0.0) $
                              , upper = 3, lower = 3 $
                              , bkspace = bsp, maxiter = 20, maxrej = 10)
       sky_resids[inslit] = bspline_valu(piximg[inslit], sset)
       IF KEYWORD_SET(CHK) THEN BEGIN
          plotx = piximg[fitpix[psort]]
          ploty = obj_min_sky[fitpix[psort]]
          rms = sqrt(djs_median(ploty^2))
          x_splot, plotx, ploty, psym1 = 3 $
                   , ymnx = [-10.0*rms, 10.0*rms] $
                   , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
;                     , xrange = [min(plotx), max(plotx)] $
       ENDIF
    ENDFOR
 ENDIF
 
 
 IF keyword_set(CHK) THEN BEGIN
    xatv, (obj_min_sky-sky_resids)*sqrt(ivar), /block
 ENDIF
 sky_model = sky_img + sky_resids

 RETURN, obj_min_sky
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
