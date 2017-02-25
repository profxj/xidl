PRO QSO_SMG_OFFSET, targname, scipath, flip = flip $
                    , hand_x = hand_x $
                    , hand_y = hand_y, hand_fwhm = hand_fwhm_out $
                    , x_qso = x_qso $
                    , scifiles = scifiles

  IF KEYWORD_SET(FLIP) THEN SIGN = -1.0d ELSE SIGN = 1.0d
  smg_file = '/Users/joe/Dropbox/SMGs/lris.1312/sXb_xID_obs.fits'
  smg = mrdfits(smg_file, 1)
  itarg = WHERE(strmatch(smg.targ, '*' + targname + '*'))
  z_qso = smg[itarg].Z
  theta = 3600.0d*djs_diff_angle(smg[itarg].QSO_RA, smg[itarg].QSO_DEC $
                                 , smg[itarg].XIDRA2, smg[itarg].XIDDEC2)
  scifiles = findfile(scipath + '*.fits.gz')
  nfiles = n_elements(scifiles)
  scihdr = headfits(scifiles[0])
  bin = long(strsplit(sxpar(scihdr, 'BINNING'), ',', /extract))
  ybin = bin[1]
  plate_scale = ybin*0.135d
  ;; This is a kludge to distinguish red and blue
  IF ybin EQ 1 THEN lam_lya = 7500.0 $ 
  ELSE IF ybin EQ 2 THEN lam_lya = 4500.0 ;;(1.0d + z_qso)*1215.67d
  

  hand_y = fltarr(nfiles)
  hand_x = fltarr(nfiles)
  hand_fwhm = fltarr(nfiles)
  x_qso = fltarr(nfiles)
  FOR ii = 0L, nfiles-1L DO BEGIN
     obj = mrdfits(scifiles[ii], 5)
     peak = max(obj.PEAKFLUX, iqso)
     nspec = n_elements(obj[iqso].WAVE_OPT)
     hand_fwhm[ii] = median(obj[iqso].FWHMFIT) 
     hand_y[ii] = interpol(findgen(nspec), obj[iqso].WAVE_OPT, lam_lya)
     x_qso[ii] = interpol(obj[iqso].XPOS, obj[iqso].YPOS, hand_y[ii])
     hand_x[ii] = x_qso[ii] - SIGN*(theta/plate_scale)
  ENDFOR
;mom1 = moment(hand_y)
;hand_y_out = replicate(mom1[0], nfiles)
  mom2 = moment(hand_fwhm)
  hand_fwhm_out = replicate(mom2[0], nfiles)
  
  print, '           FILE   X_QSO   HAND_X   HAND_Y HAND_FWHM'
  forprint, repstr(fileandpath(scifiles), '.fits.gz') $
            , x_qso, hand_x, hand_y, hand_fwhm_out $
            , textout = 1, format = 'A,F8.1,F8.1,F8.1,F6.2'
RETURN
END
