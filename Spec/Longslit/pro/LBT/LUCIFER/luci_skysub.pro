FUNCTION LUCI_SKYSUB, afiles, bfiles, tset_slits $
                      , ivar = ivar_bar, waveimg = waveimg_bar $
                      , STDTRACE = STDTRACE $
                      , sky = sky_bar $
                      , obj_pos = obj_pos, obj_neg = obj_neg $ 
                      , TELLURIC = TELLURIC1, chk = chk, wvchk = wvchk $
                      , qafile = qafile $
                      , pixflatfile = pixflatfile $
                      , illumflatfile = illumflatfile, darkfile = darkfile $
                      , SOFI = SOFI, WAVEFILE_TELL = WAVEFILE_TELL $
                      , CALIB = CALIB, SIZE_OBJMASK = SIZE_OBJMASK $
                      , PEAKTHRESH = PEAKTHRESH

IF KEYWORD_SET(TELLURIC1) THEN TELLURIC = TELLURIC1 ELSE TELLURIC = 0

dims = tset_slits[0].DIMS
nx = dims[0]
ny = dims[1]
nseq = n_elements(afiles)

FOR ii = 0L, nseq-1L DO BEGIN
   ;; Read in a and b files
   IF KEYWORD_SET(SOFI) THEN BEGIN
      sofi_proc, afiles[ii], aimg, ivar_a, hdr = hdr_a, pixflatfile = pixflatfile $
                 , illumflatfile = illumflatfile, darkfile = darkfile
      sofi_proc, bfiles[ii], bimg, ivar_b, hdr = hdr_b, pixflatfile = pixflatfile $
                 , illumflatfile = illumflatfile, darkfile = darkfile
   ENDIF ELSE BEGIN 
      luci_proc, afiles[ii], aimg, ivar_a, hdr = hdr_a, pixflatfile = pixflatfile $
                 , illumflatfile = illumflatfile, darkfile = darkfile
      luci_proc, bfiles[ii], bimg, ivar_b, hdr = hdr_b, pixflatfile = pixflatfile $
                 , illumflatfile = illumflatfile, darkfile = darkfile
   ENDELSE
   ;; Do first sky-subtraction
   mask_a = (ivar_a GT 0.0)
   mask_b = (ivar_b GT 0.0)
   var_a = (ivar_a GT 0.0)/(ivar_a + (ivar_a LE 0.0))
   var_b = (ivar_b GT 0.0)/(ivar_b + (ivar_b LE 0.0))
   a_min_b = aimg-bimg
   var_a_min_b = var_a + var_b
   mask_ab = (mask_a GT 0.0) AND (mask_b GT 0.0)
   ivar = float(mask_ab)/(var_a_min_b + (var_a_min_b EQ 0.0))
   avg_sky = (aimg + bimg)/2.0
   IF KEYWORD_SET(TELLURIC) THEN BEGIN
      ;;IF TELLURIC[0] EQ 'use_tell_wave' THEN BEGIN
      ;;   arcimg = avg_sky
      ;;   arcivar = ivar
      ;;   hdr_arc = hdr_a
      ;;ENDIF ELSE
      IF KEYWORD_SET(SOFI) THEN sofi_proc, WAVEFILE_TELL, arcimg, arcivar, hdr = hdr_arc $
                                           , pixflatfile = pixflatfile $
                                           , illumflatfile = illumflatfile, darkfile = darkfile $
      ELSE luci_proc, WAVEFILE_TELL, arcimg, arcivar, hdr = hdr_arc $
                      , pixflatfile = pixflatfile $
                      , illumflatfile = illumflatfile, darkfile = darkfile
   ENDIF ELSE BEGIN
      arcimg = avg_sky
      arcivar = ivar
      hdr_arc = hdr_a
   ENDELSE
   waveimg = luci_waveimg(arcimg, arcivar, tset_slits, hdr_arc $
                          , piximg = piximg, CHK = WVCHK, QAFILE = QAFILE $
                          , CALIB = CALIB)

   IF NOT KEYWORD_SET(TELLURIC) THEN $
      sky_resids = luci_skysub_ab(a_min_b, ivar, tset_slits, piximg $
                                  , FWHM = FWHM, CHK = CHK, SOFI = SOFI $
                                  , STDTRACE = STDTRACE, SIZE_OBJMASK = SIZE_OBJMASK $
                                  , PEAKTHRESH = PEAKTHRESH) $
   ELSE sky_resids = 0.0*a_min_b
   
   IF ii EQ 0 THEN BEGIN
      a_min_b_stk = fltarr(nx, ny, nseq)
      sky_resids_stk = fltarr(nx, ny, nseq)
      ivar_stk    = fltarr(nx, ny, nseq)
      mask_stk = fltarr(nx, ny, nseq)
      piximg_stk = fltarr(nx, ny, nseq)
      waveimg_stk = fltarr(nx, ny, nseq)
      sky_stk = fltarr(nx, ny, nseq)
      ;;hdr = hdr_a
      ;;arc_ref = luci_extract_arc(aimg, ivar_a, trace)
      ;;log_arc_ref = rebin(alog10(arc_ref > 1.0), ny*nsamp)
   ENDIF
   a_min_b_stk[*, *, ii] = a_min_b
   sky_resids_stk[*, *, ii] = sky_resids
   ivar_stk[*, *, ii] = ivar
   mask_stk[*, *, ii] = mask_ab
   piximg_stk[*, *, ii] = piximg
   waveimg_stk[*, *, ii] = waveimg
   sky_stk[*, *, ii] = avg_sky
ENDFOR
var_stk = (ivar_stk GT 0.0)/(ivar_stk + (ivar_stk LE 0.0))

diff_stk = a_min_b_stk-sky_resids_stk
IF nseq GT 1 THEN BEGIN
   ;; Avsigclip to generate a mask, i.e. CR and hot pixel reject
   diff_avs = djs_avsigclip(diff_stk, 3 $
                            , sigrej = sigrej, maxiter = maxiter $
                            , inmask = (mask_stk EQ 0) $
                            , outmask = outmask)
   weights = float(outmask EQ 0)
   w_sum = total(weights, 3)
   diff_bar = (w_sum GT 0)*total(weights*diff_stk, 3)/(w_sum + (w_sum EQ 0.0))
   var_bar = (w_sum GT 0.0)*total(weights^2*var_stk, 3)/(w_sum + $
                                                         (w_sum EQ 0.0))^2
   ivar_bar = (var_bar GT 0.0)/(var_bar + (var_bar LE 0.0))

   waveimg_bar = (w_sum GT 0)*total(weights*waveimg_stk, 3)/(w_sum + (w_sum EQ 0.0))
   sky_bar = (w_sum GT 0)*total(weights*sky_stk, 3)/(w_sum + (w_sum EQ 0.0))
ENDIF ELSE BEGIN
   diff_bar = diff_stk*mask_stk
   ivar_bar = ivar_stk*mask_stk
   waveimg_bar = waveimg_stk
   sky_bar = sky_stk
ENDELSE

;; Now do object finding on combined high S/N ratio stacked image
print, "Finding objects in stacked sky-subtracted image, first pass"
; SOFI plate scale is 0.288. Assuming 1.0" seeing FWHM = 3.5 pix
obj_pos = long_objfind(diff_bar, tset_slits = tset_slits $
                       , FWHM = FWHM $
                       , NPERSLIT = KEYWORD_SET(TELLURIC) $
                       , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                       , HAND_X = HAND_X, HAND_Y = HAND_Y $
                       , HAND_FWHM = HAND_FWHM, STDTRACE = STDTRACE $
                       , peakthresh = peakthresh)
IF KEYWORD_SET(OBJ_POS) THEN fwhm = total(obj_pos.FWHM)/n_elements(obj_pos.FWHM)
; Search for negative objects that might be left over in the avg sky
obj_neg = long_objfind(-diff_bar, tset_slits = tset_slits $
                       , FWHM = fwhm $
                       , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                       , NPERSLIT = KEYWORD_SET(TELLURIC) $
                       , HAND_X = HAND_X, HAND_Y = HAND_Y $
                       , HAND_FWHM = HAND_FWHM, STDTRACE = STDTRACE $
                       , peakthresh = peakthresh)


                       
IF NOT KEYWORD_SET(TELLURIC) THEN BEGIN 
   ;; If there is more than one sequence at the same position, refine the
   ;; sky-subtraction and object finding, but stacking the sky-subtracted frames
   IF nseq GT 1 THEN BEGIN 
      ;;OBJMASK = (OBJMASK_POS EQ 1) OR (OBJMASK_NEG EQ 1)
      ;; Update the sky-subtraction with these object locations
      FOR ii = 0L, nseq-1L DO BEGIN
         sky_resids_stk[*, *, ii] = luci_skysub_ab(a_min_b_stk[*, *, ii] $
                                                   , ivar_stk[*, *, ii] $
                                                   , tset_slits $
                                                   , piximg_stk[*, *, ii] $
                                                   , CHK = CHK $
                                                   , obj_pos = obj_pos $
                                                   , obj_neg = obj_neg $
                                                   , FWHM = FWHM $
                                                   , SOFI = SOFI $
                                                   , STDTRACE = STDTRACE $
                                                   , SIZE_OBJMASK = SIZE_OBJMASK $
                                                   , peakthresh = peakthresh)
      ENDFOR
      ;; Recombine with this modified sky-subtraction 
      diff_stk = a_min_b_stk-sky_resids_stk
      ;; Avsigclip to generate a mask, i.e. CR and hot pixel reject
      diff_avs = djs_avsigclip(diff_stk, 3 $
                               , sigrej = sigrej, maxiter = maxiter $
                               , inmask = (mask_stk EQ 0) $
                               , outmask = outmask)
      weights = float(outmask EQ 0)
      w_sum = total(weights, 3)
      diff_bar = (w_sum GT 0)*total(weights*diff_stk, 3)/(w_sum + (w_sum EQ 0.0))
      var_bar = (w_sum GT 0.0)*total(weights^2*var_stk, 3)/(w_sum + $
                                                            (w_sum EQ 0.0))^2
      ivar_bar = (var_bar GT 0.0)/(var_bar + (var_bar LE 0.0))
      
      waveimg_bar = (w_sum GT 0)*total(weights*waveimg_stk, 3)/(w_sum + (w_sum EQ 0.0))
      sky_bar = (w_sum GT 0)*total(weights*sky_stk, 3)/(w_sum + (w_sum EQ 0.0))
;; Update the object finding one last time
      print, "Finding objects in stacked sky-subtracted image, second pass"
      obj_pos = long_objfind(diff_bar, tset_slits = tset_slits $
                             , FWHM = FWHM $
                             , NPERSLIT = KEYWORD_SET(TELLURIC) $
                             , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                             , HAND_X = HAND_X, HAND_Y = HAND_Y $
                             , HAND_FWHM = HAND_FWHM, STDTRACE = STDTRACE $
                             , peakthresh = peakthresh)

      IF KEYWORD_SET(OBJ_POS) THEN fwhm = total(obj_pos.FWHM)/n_elements(obj_pos.FWHM)
      obj_neg = long_objfind(-diff_bar, tset_slits = tset_slits $
                             , FWHM = fwhm $
                             , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                             , NPERSLIT = KEYWORD_SET(TELLURIC) $
                             , HAND_X = HAND_X, HAND_Y = HAND_Y $
                             , HAND_FWHM = HAND_FWHM, STDTRACE = STDTRACE $
                             , peakthresh = peakthresh)
   ENDIF
ENDIF
   
IF KEYWORD_SET(CHK) THEN xatv, diff_bar*sqrt(ivar_bar), wv = waveimg_bar, block = chk

RETURN, diff_bar
END
