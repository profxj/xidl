FUNCTION GNIRS_SKYSUB, filenames, flatfile, acqfiles, OBJMASK = OBJMASK1 $
                       , IVAR_ABBA = IVAR_ABBA $
                       , OBJ_POS = OBJ_POS, OBJ_NEG = OBJ_NEG $
                       , WAVEIMG = WAVEIMG, SKY_MODEL = SKY_MODEL $
                       , TELLURIC = TELLURIC1, ZAP = ZAP1, CHK = CHK $
                       , VERBOSE = VERBOSE, HDR = HDR, targdir = targdir $
                       , MINFWHM = MINFWHM, MAXFWHM = MAXFWHM $
                       , SKY_RESIDS = SKY_RESIDS, WVCHK = WVCHK

;, IVAR_POS = IVAR_POS, 
;, IVAR_NEG = IVAR_NEG, 

;;  chk = 1
plate_scale = 0.15D
IF NOT KEYWORD_SET(BSP) THEN BSP = 0.8
IF NOT KEYWORD_SET(MAXFWHM) THEN MAXFWHM = 1.3D/plate_scale
;; Minimum FWHM is 1 pixel
IF NOT KEYWORD_SET(MINFWHM) THEN MINFWHM = 1.0D
IF KEYWORD_SET(ZAP1) THEN ZAP = ZAP1 $
ELSE ZAP = 0
IF KEYWORD_SET(OBJMASK1) THEN OBJMASK = OBJMASK1 $
ELSE OBJMASK = 0
IF KEYWORD_SET(TELLURIC1) THEN BEGIN
    TELLURIC = TELLURIC1
    OBJMASK = 0
ENDIF ELSE TELLURIC = 0

IF (N_params() LT 1) THEN BEGIN
    print, 'Syntax: gnirs_skysub,filenames,flatfile,[OBJFIND=, ZAP=, CHK=,VERBOSE=,]'
    return, 0
endif

flat = mrdfits(flatfile, 0)
tset_slits = mrdfits(flatfile, 1)
dimt = size(tset_slits.coeff, /dimen)
norders = dimt[1]
dims = size(flat, /dimen)
nx = dims[0]
ny = dims[1]
order_vec = [3, 4, 5, 6, 7, 8]

;; order boundaries
traceset2xy, tset_slits[0], yy1, xx1
traceset2xy, tset_slits[1], yy2, xx2
ordermask = long_slits2mask(tset_slits)
;      enumerate slits by order
ordermask[WHERE(ordermask GT 0)] = ordermask[WHERE(ordermask GT 0)] + 2L
ximg = long_slits2x(tset_slits, edgmask = edgmask)
;      Using ABBA sky subtraction.  Must use this if no darks were taken...
gnirs_abba_proc, filenames, flat, tset_slits, acqfiles, ZAP = ZAP $
                 , abba = abba, ivar_abba = ivar_abba $
                 , piximg = piximg, waveimg = waveimg, sky_abba = sky_abba $
                 , TELLURIC = TELLURIC, hdr = hdr, targdir = targdir $
                 , WVCHK = WVCHK
;; , ivar_pos = ivar_pos, ivar_neg = ivar_neg $

;   Second-pass bspline sky subtraction

IF KEYWORD_SET(TELLURIC) THEN BEGIN
;   Skip second-pass for standard stars and bright objects.
    print, "    Bright calibration source; skipping 2nd pass..."
    sky_model = sky_abba 
    sky_resids = 0.0*sky_abba
ENDIF ELSE BEGIN
;      Otherwise do 2nd pass
    print, "  Second-pass sky subtraction..."
    
    y_img = replicate(1.0, nx)#findgen(ny)
    sky_resids = abba*0.0
    FOR iorder = 0, norders-1 DO BEGIN
        inorder = where(ordermask EQ order_vec[iorder], nord)
        fitpix = where(ordermask EQ order_vec[iorder] AND $
                       finite(abba)  AND $
                       edgmask EQ 0 AND $
                       abs(abba) LE 5.0d4 AND $
                       ivar_abba GT 0.0, npix)
;                       ximg GT 0.1  AND $
;                       ximg LT 0.9  AND $
        psort = sort(piximg[fitpix])
        sset = bspline_iterfit(piximg[fitpix[psort]], abba[fitpix[psort]] $
                               , invvar = (ivar_abba[fitpix[psort]] GT 0.0) $
                               , upper = 3, lower = 3 $
                               , bkspace = bsp, maxiter = 20, maxrej = 10 $
                               , /silent)
        sky_resids[inorder] = bspline_valu(piximg[inorder] $
                                           , sset)
        sky_model = sky_abba + sky_resids
        IF KEYWORD_SET(CHK) THEN BEGIN
            plotx = piximg[fitpix[psort]]
            ploty = abba[fitpix[psort]]
            rms = sqrt(djs_median(ploty^2))
            x_splot, plotx, ploty, psym1 = 3 $
                     , ymnx = [-10.0*rms, 10.0*rms] $
                     , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
         ENDIF
     ENDFOR
ENDELSE

print, "Finding objects in sky-subtracted image"

FWHM = 5.0
ncoeff = 5
obj_pos = long_objfind(abba-sky_resids, tset_slits = tset_slits, nperslit = 1 $
                       , FWHM = FWHM, ncoeff = ncoeff)
flag_ovlp_pos = 0
IF size(obj_pos, /type) EQ 8 THEN BEGIN 
   ipos_bad = WHERE(total(obj_pos.XPOS, 1) EQ 0, npos_bad)
   FOR iobj = 0L, n_elements(obj_pos)-1L DO BEGIN
      ipos_ovlp = WHERE(obj_pos[iobj].XPOS LT xx1[*, obj_pos[iobj].SLITID-1L] OR $
                        obj_pos[iobj].XPOS GT xx2[*, obj_pos[iobj].SLITID-1L] $
                        , npos_ovlp)
      IF npos_ovlp GT 0 THEN flag_ovlp_pos = 1
   ENDFOR
ENDIF ELSE npos_bad = 6
   

IF n_elements(obj_pos) NE 6 OR npos_bad GT 0 OR flag_ovlp_pos THEN BEGIN
    print, 'Failed to find 6 good positive objects ... reiterating with smaller FWHM'
    obj_pos = long_objfind(abba-sky_resids, tset_slits = tset_slits $
                           , nperslit = 1, FWHM = FWHM/2.0, ncoeff = ncoeff)
ENDIF
obj_pos.SLITID = obj_pos.SLITID+2L

fwhm = total(obj_pos.FWHM)/n_elements(obj_pos.FWHM)
obj_neg = long_objfind(sky_resids-abba, tset_slits = tset_slits, nperslit = 1 $
                       , FWHM = fwhm, ncoeff = ncoeff)

flag_ovlp_neg = 0
IF size(obj_neg, /type) EQ 8 THEN BEGIN 
   ineg_bad = WHERE(total(obj_neg.XPOS, 1) EQ 0, nneg_bad)
   FOR iobj = 0L, n_elements(obj_neg)-1L DO BEGIN
      ineg_ovlp = WHERE(obj_neg[iobj].XPOS LT xx1[*, obj_neg[iobj].SLITID-1L] OR $
                        obj_neg[iobj].XPOS GT xx2[*, obj_neg[iobj].SLITID-1L] $
                        , nneg_ovlp)
      IF nneg_ovlp GT 0 THEN flag_ovlp_neg = 1
   ENDFOR
ENDIF ELSE nneg_bad = 6

IF n_elements(obj_neg) NE 6 OR nneg_bad GT 0 OR flag_ovlp_neg THEN BEGIN
   print, 'Failed to find 6 negative objects ... reiterating with smaller FWHM'
   obj_neg = long_objfind(sky_resids-abba, tset_slits = tset_slits $
                          , nperslit = 1, FWHM = fwhm/2.0, ncoeff = ncoeff)
ENDIF
obj_neg.SLITID = obj_neg.SLITID+2L

IF KEYWORD_SET(OBJMASK) THEN BEGIN
    print, "Doing third pass sky subtraction with object masking"
    pos_fwhm = (obj_pos.FWHM >  MINFWHM) <  MAXFWHM
    neg_fwhm = (obj_neg.FWHM >  MINFWHM) <  MAXFWHM
    left_pos = obj_pos.xpos -  $
      0.5*obj_pos.FWHM ## replicate(1.0D, ny)
    right_pos = obj_pos.xpos + $
      0.5*obj_pos.FWHM ## replicate(1.0D, ny)
    left_neg = obj_neg.xpos -  $
      0.5*obj_neg.FWHM ## replicate(1.0D, ny)
    right_neg = obj_neg.xpos + $
      0.5*obj_neg.FWHM ## replicate(1.0D, ny)
;                   xmap = findgen(nx) # replicate(1.0, ny)
    skymask = long(0.0*abba)
    skymask[WHERE(ordermask)] = 1L
    FOR iorder = 0L, norders-1L DO BEGIN
        FOR j = 0L, ny-1L DO BEGIN 
            xmin_pos = floor(left_pos[j, iorder]) >  0
            xmax_pos = ceil(right_pos[j, iorder]) < (nx-1L)
            xmin_neg = floor(left_neg[j, iorder]) >  0
            xmax_neg = ceil(right_neg[j, iorder]) < (nx-1L)
            skymask[xmin_pos:xmax_pos, j] = 0L
            skymask[xmin_neg:xmax_neg, j] = 0L
        ENDFOR
    ENDFOR
    sky_resids2 = abba*0.0
    FOR iorder = 0, norders-1 DO BEGIN
        buffer = 4
        inorder = where(ordermask EQ order_vec[iorder], nord)
        fitpix = where(ordermask EQ order_vec[iorder] AND $
                       finite(abba)     AND $
                       ivar_abba GT 0.0 AND $
                       edgmask EQ 0    AND $
                       abs(abba) LE 5.0d4 AND $
                       skymask EQ 1, npix)
;                       ximg GT 0.1     AND $
;                       ximg LT 0.9     AND $
;                       y_img LT 1005    AND $
;        med_rms = sqrt(djs_median(abba[fitpix]^2))
;        fitpix = fitpix[WHERE(abba[fitpix] LE 20.0*med_rms AND  $
;                              med_rms GE -20.0*med_rms)]
        psort = sort(piximg[fitpix])
        sset = bspline_iterfit(piximg[fitpix[psort]], abba[fitpix[psort]] $
                               , invvar = (ivar_abba[fitpix[psort]] GT 0.0) $
                               , upper = 3, lower = 3 $
                               , bkspace = bsp, maxiter = 20, maxrej = 10)
        sky_resids2[inorder] = bspline_valu(piximg[inorder] $
                                           , sset)
        IF KEYWORD_SET(CHK) THEN BEGIN
            plotx = piximg[fitpix[psort]]
            ploty = abba[fitpix[psort]]
            rms = sqrt(djs_median(ploty^2))
            x_splot, plotx, ploty, psym1 = 3 $
                     , ymnx = [-5.0*rms, 5.0*rms] $
                     , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
        ENDIF
    ENDFOR
    sky_resids = sky_resids2
    sky_model = sky_abba + sky_resids

 ENDIF

IF NOT KEYWORD_SET(TELLURIC) THEN BEGIN 
   ;; One more CR rejection to mask any pixels which we missed. This
   ;; causes problems with telluric since the telluric absorption
   ;; features look like cosmics
   FWHM = 3.0                  
   sigma_psf = FWHM/2.35482D
   psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
   cr_sharp_pos = psf_reject_cr(abba - sky_resids, ivar_abba, psfvals $
                                , satmask = (abba-sky_resids LT -10))
   cr_sharp_neg = psf_reject_cr(sky_resids-abba, ivar_abba, psfvals $
                                , satmask = (sky_resids-abba LT -10))
   finalmask = (cr_sharp_pos EQ 0) AND (cr_sharp_neg EQ 0)
   ivar_abba = finalmask*ivar_abba
ENDIF

IF keyword_set(CHK) THEN $
   xatv, (abba-sky_resids)*sqrt(ivar_abba), wv = waveimg, /block
RETURN, abba
END
