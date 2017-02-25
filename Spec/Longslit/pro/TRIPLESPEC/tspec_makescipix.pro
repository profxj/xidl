;
; INPUTS: arcimg: a processed 2D arc frame or skyline spectrum
;         arc2:   (optional): if arcimg is a skyline spectrum,
;                 user may also input an additional ThAr arc spectrum
;                 to fill in orders that are sparse in sky lines
;
; RETURNS: 2D pixel image for sky subtraction
;          (optional) pixset is the tset used to generate the fit.
;
;
FUNCTION tspec_makescipix, arcimg, tset_slits, chk = chk $
                           , verbose = verbose, std = std $
                           , bright = bright, pixset = pixset $
                           , _EXTRA = keys
  
;  IF NOT KEYWORD_SET(PKWKDTH) THEN PKWDTH=5.0D
;  IF NOT KEYWORD_SET(TOLER) THEN TOLER=2.0D
;  IF NOT KEYWORD_SET(SIG_THRESH) THEN SIG_THRESH=15.0D ; 15 default
;  IF NOT KEYWORD_SET(NSIG) THEN NSIG=15.0d ; 15 default
;  IF NOT KEYWORD_SET(FWHM) THEN FWHM=3.0D
;  IF NOT KEYWORD_SET(BOX_RADIUS) THEN BOX_RADIUS=5.0D
;  IF NOT KEYWORD_SET(med_err) THEN med_err = 0.16D
  
  ;;omask_trim = tspec_ordermask(tset_slits, /fsr)
  omask = tspec_ordermask(tset_slits)
  omask_trim = omask
  imask = (arcimg GT -20.0 AND arcimg LT 1d5 AND omask_trim GT 0)

  ; Case where only ONE frame (e.g. sky) used to determine tilts.
;  if (NOT keyword_set(ARC2)) then begin
  chk = 1
  wset = tspec_wavepix(imask*arcimg, tset_slits, fwhm = fwhm $    
                       , sig_thresh = sig_thresh $
                       , nsig = NSIG $
                       , pkwdth = PKWDTH $
                       , TOLER = TOLER $
                       , med_err = med_err $
                       , BOX_RADIUS = box_radius $
                       , CHK = chk $
                       , verbose = verbose $
                       , _EXTRA = keys ) 
;                              , THAR = thar $
  ;;piximg_arc = long_wpix2image(wset_arc, tset_slits)
  
;;  endif else begin     ; Case where multiple frames used (default in fire_pipe)
     
;     restore, getenv("FIRE_DIR")+'/Calib/fire_piximg_wset.idl'
;     piximg_in = long_wpix2image(wset, tset_slits)
     
;     imask = (arcimg GT -20.0 AND arcimg LT 1d5)
;     wset_arc = fire_wavepix(imask*(arcimg+arc2), tset_slits, fwhm=fwhm $    
;                             , sig_thresh = sig_thresh $
;                             , nsig=NSIG $
;                             , pkwdth = PKWDTH $
;                             , TOLER = TOLER $
;                             , med_err = med_err $
;                             , CHK = chk $
;                             , BOX_RADIUS=box_radius $
;                             , verbose = verbose $
;;                             , piximg_in = piximg_in $
;                             , THAR=thar $
;                             , ONLY_SLITS = [1,2,3,4,5,6,7,8,9,10,11,12,13,$
;                                             14,15,16,17,18,19,20,21] $
;                             , BAD_SLITS=badslits, _EXTRA = keys ) 

 ;    wset_arc = fire_wavepix(imask*arcimg, tset_slits $
 ;                            , CHK = chk $
 ;                            , THAR=thar $
 ;                            , ONLY_SLITS = [1,2,3,4,5,6,7,8,9,10,11,12,13,$
 ;                                            14,15,16,17,18,19,20,21] $
 ;                            , BAD_SLITS=badslits, _EXTRA = keys ) 
 ; endelse
  
     ;;wset = wset_arc

  ; Only fall back on the ThAr spectrum for orders with bad pixel fits.
  ; This typically happens in the bluest few orders.
;  if (keyword_set(arc2) AND n_elements(badslits) GT 1 AND NOT keyword_set(THAR)) then begin
;     sciimg = arc2
;     imask = (sciimg GT -20.0 AND sciimg LT 1d5)
;     
;     if (n_elements(slits) GT 0) then begin
;        slits = [slits, (badslits+1)]
;     endif else begin
;        slits = [(badslits+1)]
;     endelse 

;     stop
     ; Don't mask out the ends of the orders for 2nd pass.
;     wset_sky = fire_wavepix(imask*sciimg, tset_slits, FWHM = FWHM $
;                             , pkwdth = pkwdth, toler = toler $
;                             , sig_thresh=sig_thresh $
;                             , nsig = nsig $
;;                             , ONLY_SLITS = slits $
;                             , BOX_RADIUS=box_radius $
;;;                              , piximg_in = piximg_arc $
;                             , med_err = med_err $
;                             , CHK=chk $
;                             , THAR=thar $
;                             , BAD_SLITS=badslits2, _EXTRA = keys )
;     slitind = slits-1
;     wset[slitind] = wset_sky[slitind]
;  endif

;  wset_new = fire_wset_clean(wset, imask*arcimg, tset_slits)

  piximg = long_wpix2image(wset, tset_slits)

  RETURN,piximg
  
end
