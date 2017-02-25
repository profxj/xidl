FUNCTION mage_makescipix, arcimg, sciimg, tset_slits, chk=chk,std=std, bright=bright
  
  IF NOT KEYWORD_SET(PKWKDTH) THEN PKWDTH=3.0D
  IF NOT KEYWORD_SET(TOLER) THEN TOLER=2.0D
  IF NOT KEYWORD_SET(SIG_THRESH) THEN SIG_THRESH=5.0D
  IF NOT KEYWORD_SET(NSIG) THEN NSIG=5.0d
  IF NOT KEYWORD_SET(FWHM) THEN FWHM=3.0D
  IF NOT KEYWORD_SET(BOX_RADIUS) THEN BOX_RADIUS=3.0D

  wset_arc = long_wavepix(arcimg, tset_slits, fwhm=fwhm $
                          , box_radius = box_radius $
                          , sig_thresh = sig_thresh $
                          , nsig=NSIG $
                          , pkwdth = PKWDTH $
                          , TOLER = TOLER $
                          , CHK = CHK) 
  piximg_arc = long_wpix2image(wset_arc, tset_slits)
  IF KEYWORD_SET(STD) OR KEYWORD_SET(BRIGHT) THEN piximg=piximg_arc $
  ELSE BEGIN
     imask = (sciimg GT -20.0 AND sciimg LT 1d5)
     wset_sky = long_wavepix(imask*sciimg, tset_slits, FWHM = FWHM $
                             , pkwdth = pkwdth, toler = toler $
                             , sig_thresh=sig_thresh $
                             , nsig = nsig $
                             , ISLIT = [12,13,14,15] $
                             , piximg_in = piximg_arc $
                             , CHK=CHK)
     wset = wset_arc
     wset[11:*] = wset_sky[11:*]
     piximg = long_wpix2image(wset, tset_slits)
  ENDELSE
  RETURN,piximg

end
