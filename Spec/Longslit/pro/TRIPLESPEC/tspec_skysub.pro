FUNCTION TSPEC_SKYSUB, filenames, orderfile, pixflatfile, illumflatfile $
                       , tset_slits $
                       , OBJMASK = OBJMASK1 $
                       , IVAR_AB = IVAR_AB $
                       , OBJ_POS = OBJ_POS, OBJ_NEG = OBJ_NEG $
                       , WAVEIMG = WAVEIMG, SKY_MODEL = SKY_MODEL $
                       , SCIHDR = SCIHDR $
                       , TELLURIC = TELLURIC1, ZAP = ZAP1, CHK = CHK $
                       , VERBOSE = VERBOSE, HDR = HDR, targdir = targdir $
                       , MINFWHM = MINFWHM, MAXFWHM = MAXFWHM $
                       , SKY_RESIDS = SKY_RESIDS, WVCHK = WVCHK $
                       , EX_BOTH=ex_both

IF NOT KEYWORD_SET(SKYBUFFER) THEN SKYBUFFER = 30L
plate_scale = 0.234D            ; TSPEC plate scale
IF NOT KEYWORD_SET(BSP) THEN BSP = 0.5
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
    print, 'Syntax: tspec_skysub,filenames,flatfile,[OBJFIND=, ZAP=, CHK=,VERBOSE=,]'
    return, 0
endif

dimt = size(tset_slits.coeff, /dimen)
ordermask = tspec_ordermask(tset_slits, order_vec = order_vec)
norders = dimt[1]
dims = size(ordermask, /dimen)
nx = dims[0]
ny = dims[1]
;; order boundaries
traceset2xy, tset_slits[0], yy1, xx1
traceset2xy, tset_slits[1], yy2, xx2
ximg = long_slits2x(tset_slits, edgmask = edgmask)
;; Using AB sky subtraction.  
TSPEC_DIFF_PROC, filenames, pixflatfile, illumflatfile $
                 , tset_slits, ZAP = ZAP $
                 , AB = AB, ivar_AB = ivar_AB, piximg = piximg $
                 , targdir = targdir, sky_AB = sky_AB $
                 , waveimg = waveimg, TELLURIC = TELLURIC, hdr = scihdr $
                 , WVCHK = WVCHK
;   Second-pass bspline sky subtraction
IF KEYWORD_SET(TELLURIC) THEN BEGIN
;   Skip second-pass for standard stars and bright objects.
    print, "    Bright calibration source; skipping 2nd pass..."
    sky_model = sky_AB 
    sky_resids = 0.0*sky_AB
ENDIF ELSE BEGIN
;      Otherwise do 2nd pass
   print, "  Second-pass sky subtraction..."
   y_img = replicate(1.0, nx)#findgen(ny)
   sky_resids = AB*0.0
   FOR iorder = 0, norders-1 DO BEGIN
      inorder = where(ordermask EQ order_vec[iorder], nord)
      fitpix = where(ordermask EQ order_vec[iorder] AND $
                     finite(AB)  AND $
                     edgmask EQ 0 AND $
                     abs(AB) LE 5.0d4 AND $
                     ivar_AB GT 0.0, npix)
;                       ximg GT 0.1  AND $
;                       ximg LT 0.9  AND $
      psort = sort(piximg[fitpix])
      sset = bspline_iterfit(piximg[fitpix[psort]], AB[fitpix[psort]] $
                             , invvar = (ivar_AB[fitpix[psort]] GT 0.0) $
                             , upper = 3, lower = 3 $
                             , bkspace = bsp, maxiter = 20, maxrej = 10 $
                             , /silent)
      sky_resids[inorder] = bspline_valu(piximg[inorder], sset)
      sky_model = sky_AB + sky_resids
      IF KEYWORD_SET(CHK) THEN BEGIN
         plotx = piximg[fitpix[psort]]
         ploty = AB[fitpix[psort]]
         rms = sqrt(djs_median(ploty^2))
         ;x_splot, plotx, ploty, psym1 = 3 $
         ;         , ymnx = [-10.0*rms, 10.0*rms] $
         ;         , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
      ENDIF
   ENDFOR
ENDELSE

print, "Finding objects in sky-subtracted image"
;; Run first time 
FWHM = 5.0
;stop
obj_neg = tspec_findobj(sky_resids-AB, ivar_AB, waveimg, tset_slits $
                        , filstd = filstd, CHK = chk $
                        , APERV = aperv, APERS = apers $
                        , NFIND = NFIND, PEAKTHRESH = PEAKTHRESH $
                        , FWHM = FWHM $
                        , ABSTHRESH = ABSTHRESH, MIN_SN = MIN_SN)
  if size(obj_neg,/type) EQ 2 then begin
    print, 'tspec_skysub:  No object.  Punting'
    return, -1
  endif

;obj_pos = long_objfind(AB-sky_resids, tset_slits = tset_slits, nperslit = 1 $
;                       , FWHM = FWHM)
;ipos_bad = WHERE(total(obj_pos.XPOS, 1) EQ 0, npos_bad)
;flag_ovlp_pos = 0
;FOR iobj = 0L, n_elements(obj_pos)-1L DO BEGIN
;   ipos_ovlp = WHERE(obj_pos[iobj].XPOS LT xx1[*, obj_pos[iobj].SLITID-1L] OR $
;                     obj_pos[iobj].XPOS GT xx2[*, obj_pos[iobj].SLITID-1L] $
;                     , npos_ovlp)
;    IF npos_ovlp GT 0 THEN flag_ovlp_pos = 1
;ENDFOR

;IF n_elements(obj_pos) NE NORDERS OR npos_bad GT 0 OR flag_ovlp_pos THEN BEGIN
;   print, 'Failed to find ' + string(norders) + 'good positive objects ... reit;erating with smaller FWHM'
;   obj_pos = long_objfind(AB-sky_resids, tset_slits = tset_slits $
;                          , nperslit = 1, FWHM = FWHM/2.0)
;ENDIF
obj_neg.SLITID = order_vec[obj_neg.slitid-1L]
fwhm = total(obj_neg.FWHM)/n_elements(obj_neg.FWHM)
obj_pos = tspec_findobj(AB-sky_resids, ivar_AB, waveimg, tset_slits $
                        , filstd = filstd, CHK = chk $
                        , APERV = aperv, APERS = apers $
                        , NFIND = NFIND, PEAKTHRESH = PEAKTHRESH $
                        , FWHM = FWHM $
                        , ABSTHRESH = ABSTHRESH, MIN_SN = MIN_SN)
if size(obj_pos,/type) EQ 2 then begin
   print,'tspec_skysub: No object found positive.  Punting!'
   RETURN, -1
endif
obj_pos.SLITID = order_vec[obj_pos.SLITID-1L]
;obj_neg = long_objfind(sky_resids-AB, tset_slits = tset_slits, nperslit = 1 $
;                       , FWHM = fwhm)
;ineg_bad = WHERE(total(obj_neg.XPOS, 1) EQ 0, nneg_bad)
;flag_ovlp_neg = 0
;FOR iobj = 0L, n_elements(obj_neg)-1L DO BEGIN
;   ineg_ovlp = WHERE(obj_neg[iobj].XPOS LT xx1[*, obj_neg[iobj].SLITID-1L] OR $
;                     obj_neg[iobj].XPOS GT xx2[*, obj_neg[iobj].SLITID-1L] $
;                     , nneg_ovlp)
;   IF nneg_ovlp GT 0 THEN flag_ovlp_neg = 1
;ENDFOR
;'IF n_elements(obj_neg) NE NORDERS OR nneg_bad GT 0 OR flag_ovlp_neg THEN BEGIN
;'   print, 'Failed to find ' + string(norders) + 'good positive objects ... rei;t'erating with smaller FWHM'
;'   obj_neg = long_objfind(sky_resids-AB, tset_slits = tset_slits $
;                          , nperslit = 1, FWHM = fwhm/2.0)
;ENDIF
;
;CHK=1
;stop

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

    ;for tt=0,1 do begin ;; Both objects
    for tt=0,keyword_set(EX_BOTH) do begin ;; Now only doing the positive object
       case tt of 
          0: begin
             left_sky  = obj_pos.xpos - SKYBUFFER 
             right_sky = obj_pos.xpos + SKYBUFFER 
          end
          1: begin
             left_sky  = obj_neg.xpos - SKYBUFFER 
             right_sky = obj_neg.xpos + SKYBUFFER 
          end
          else: stop
       endcase
       bufmask  = long(0.0*AB)
       objmask1 = long(0.0*AB)
;;    skymask[WHERE(ordermask)] = 1L
       FOR iorder = 0L, norders-1L DO BEGIN
          FOR j = 0L, ny-1L DO BEGIN 
             ;; Only use a region within SKYBUFFER of pos trace for sky-subt
             xmin_sky = floor(left_sky[j, iorder]) >  0
             xmax_sky = ceil(right_sky[j, iorder]) < (nx-1L)
             bufmask[xmin_sky:xmax_sky, j] = 1L
             ;; Mask regions near positive and negative traces
             xmin_pos = floor(left_pos[j, iorder]) >  0
             xmax_pos = ceil(right_pos[j, iorder]) < (nx-1L)
             xmin_neg = floor(left_neg[j, iorder]) >  0
             xmax_neg = ceil(right_neg[j, iorder]) < (nx-1L)
             objmask1[xmin_pos:xmax_pos, j] = 1L
             objmask1[xmin_neg:xmax_neg, j] = 1L
          ENDFOR
       ENDFOR
                                ;skymask = (ordermask GT 0) AND (objmask1 EQ 0)
       skymask = (ordermask GT 0) AND bufmask AND (objmask1 EQ 0)
       sky_resids2 = AB*0.0
       FOR iorder = 0, norders-1 DO BEGIN
          inorder = where(ordermask EQ order_vec[iorder], nord)
          fitpix = where(ordermask EQ order_vec[iorder] AND $
                         finite(AB)     AND $
                         ivar_AB GT 0.0 AND $
                         edgmask EQ 0    AND $
                         abs(AB) LE 5.0d4 AND $
                         skymask EQ 1, npix)
;                       ximg GT 0.1     AND $
;                       ximg LT 0.9     AND $
;                       y_img LT 1005    AND $
;        med_rms = sqrt(djs_median(abba[fitpix]^2))
;        fitpix = fitpix[WHERE(abba[fitpix] LE 20.0*med_rms AND  $
;                              med_rms GE -20.0*med_rms)]
          psort = sort(piximg[fitpix])
          sset = bspline_iterfit(piximg[fitpix[psort]], AB[fitpix[psort]] $
                                 , invvar = (ivar_AB[fitpix[psort]] GT 0.0) $
                                 , upper = 3, lower = 3 $
                                 , bkspace = bsp, maxiter = 20, maxrej = 10)
          sky_resids2[inorder] = bspline_valu(piximg[inorder] $
                                           , sset)
          IF KEYWORD_SET(CHK) THEN BEGIN
             plotx = piximg[fitpix[psort]]
             ploty = AB[fitpix[psort]]
             rms = sqrt(djs_median(ploty^2))
                                ;x_splot, plotx, ploty, psym1 = 3 $
                                ;        , ymnx = [-5.0*rms, 5.0*rms] $
                                ;         , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
                                ;stop
          ENDIF
       ENDFOR
       case tt of
          0: begin ;; pos obj
             sky_resids = sky_resids2
          end
          1: begin ;; replace neg obj
             gds = where(bufmask)
             sky_resids[gds] = sky_resids2[gds]
          end
       endcase
       sky_model = sky_AB + sky_resids
    endfor
 ENDIF
;; Update trace with improved sky-subtraction
fwhm = total(obj_pos.FWHM)/n_elements(obj_pos.FWHM)
obj_pos = tspec_findobj(AB-sky_resids, ivar_AB, waveimg, tset_slits $
                        , filstd = filstd, CHK = chk $
                        , APERV = aperv, APERS = apers $
                        , NFIND = NFIND, PEAKTHRESH = PEAKTHRESH $
                        , FWHM = FWHM $
                        , ABSTHRESH = ABSTHRESH, MIN_SN = MIN_SN)
if size(obj_pos,/type) EQ 2 then begin
   print,'tspec_skysub: No object found positive.  Punting!'
   RETURN, -1
endif
obj_pos.SLITID = order_vec[obj_pos.SLITID-1L]


;; One more CR rejection to mask any pixels which we missed
FWHM = 3.0                  
sigma_psf = FWHM/2.35482D
psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
cr_sharp_pos = psf_reject_cr(AB - sky_resids, ivar_AB, psfvals $
                             , satmask = (AB-sky_resids LT -10))
cr_sharp_neg = psf_reject_cr(sky_resids-AB, ivar_AB, psfvals $
                             , satmask = (sky_resids-AB LT -10))
finalmask = (cr_sharp_pos EQ 0) AND (cr_sharp_neg EQ 0)
ivar_AB = finalmask*ivar_AB

IF keyword_set(CHK) THEN $
   xatv, (AB-sky_resids)*sqrt(ivar_AB), wv = waveimg, /block

RETURN, AB
END
