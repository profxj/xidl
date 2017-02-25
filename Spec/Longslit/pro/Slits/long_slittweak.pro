FUNCTION long_slittweak, flat, tset_slits_in $
                         , CHK = CHK, GMOS = GMOS

;  IF NOT KEYWORD_SET(MINSLIT) THEN MINSLIT = 4.0
;  IF NOT KEYWORD_SET(MAXSLIT) THEN MAXSLIT = 8.0

  NORM_TOL = 0.85d
  tset_slits = tset_slits_in
  traceset2xy, tset_slits[0], rows, left_edge
  traceset2xy, tset_slits[1], rows, right_edge
  dslit = tset_slits[1].COEFF[0, *]-tset_slits[0].COEFF[0, *]
  
  slitcen = (left_edge + right_edge)/2.0d
  flat_spec = extract_boxcar(flat, slitcen, rows, radius = 2.0)/4.0
  dim = size(flat, /dim)
  nx = dim[0]
  ny = dim[1]
  xarr = findgen(nx) # replicate(1.0, ny)
  yarr = findgen(ny)## replicate(1.0, nx)
  dim = size(tset_slits[0].coeff, /dim)
  IF n_elements(dim) EQ 1 THEN BEGIN
     ncoeff = dim[0]
     nslit = 1
  ENDIF ELSE BEGIN
     ncoeff = dim[0]
     nslit = dim[1]
  ENDELSE
  flat_nrm = 0.0*flat
  shft = 2.0
  spatsamp = 0.05
  step = 0.01
  yvec = lindgen(ny)

  tset_slits_old = tset_slits
  ;; This is the default for GMOS N&S box slits
  IF KEYWORD_SET(GMOS) THEN $
     igood = where(dslit GT 2.0 AND dslit LE 30.0, ngood) $
  ELSE BEGIN
     igood = lindgen(nslit)
     ngood = nslit
  ENDELSE

  FOR ii = 0, ngood-1L DO BEGIN
     islit = igood[ii]
     flat_now1 = djs_median(flat_spec[*, islit], width = 100 $
                           , boundary = 'reflect')
     max_now = max(flat_now1)
     keep = where(flat_now1 GT 0.10*max_now)
     YSLIT_MIN = min(yvec[keep])
     YSLIT_MAX = max(yvec[keep])
     slitmask_left  = long_slits2mask(tset_slits_old, xshift = -abs(shft) $
                                      , ONLYSLITS = islit+1L)
     slitmask_right = long_slits2mask(tset_slits_old, xshift = abs(shft) $
                                      , ONLYSLITS = islit+1L)
     ipix = WHERE(slitmask_left EQ (islit + 1L) OR $
                  slitmask_right EQ (islit + 1L), npix)
     flat_now = djs_median(flat_spec[*, islit], width = 5 $
                           , boundary = 'reflect')
     flat_nrm[ipix] = flat[ipix]/interpol(flat_now, findgen(ny) $
                                          , yarr[ipix])
     ipix = WHERE((slitmask_left EQ (islit + 1L) OR $
                   slitmask_right EQ (islit + 1L)) AND $
                  flat_nrm LE 1.7 AND flat_nrm GE 0.0 AND $
                  yarr GE YSLIT_MIN AND yarr LE YSLIT_MAX, npix)
     frac = (xarr[ipix] - interpol(left_edge[*, islit], rows[*, islit], yarr[ipix]))/dslit[islit]
     psort = sort(frac)
     med_spat_width = round(double(npix)*double(spatsamp))
;    normspat = smooth(flat_nrm[ipix[psort]], med_spat_width, /edge_trunc)
      normspat = djs_median(flat_nrm[ipix[psort]] <  1.7d $
                            , width = med_spat_width, bound = 'reflect')
      ;; march out from middle to find left edge
      maxnorm = max(normspat[where(frac[psort] GE 0.4 AND frac[psort] LE 0.6D)])
      tweak_left = 0L
      FOR xleft = 0.5D, min(frac[psort]), -step DO BEGIN
         ynow = interpol(normspat, frac[psort], xleft)
         IF ynow LE NORM_TOL*maxnorm THEN BEGIN 
            tset_slits[0].COEFF[0, islit] = tset_slits[0].COEFF[0, islit] + $
                                            xleft*dslit[islit]
            tweak_left = 1L
            BREAK
         ENDIF
      ENDFOR
      tweak_right = 0L
      FOR xright = 0.5D, max(frac[psort]), step DO BEGIN
          ynow = interpol(normspat, frac[psort], xright)
          IF ynow LE NORM_TOL*maxnorm THEN BEGIN
             tset_slits[1].COEFF[0, islit] = tset_slits[1].COEFF[0, islit] - $
                                             (1.0d - xright)*dslit[islit]
             tweak_right = 1L
             BREAK
          ENDIF
       ENDFOR
      IF KEYWORD_SET(CHK) THEN BEGIN
         plot, frac[psort], flat_nrm[ipix[psort]], psym = 3 $
               , xr = [min(frac[psort]), max(frac[psort])], yr = [0.0, 2.0] $
               , title = 'slit # ' + strcompress(string(islit+1L), /rem) $
               , background = djs_icolor('black'), color = djs_icolor('white')
         oplot, frac[psort], normspat, col = djs_icolor('red'), psym = 4
         IF KEYWORD_SET(tweak_left) THEN $
            oplot, [xleft, xleft], [0.0, 2.0], linestyle = 2, col = djs_icolor('green')
         IF KEYWORd_SET(tweak_right) THEN $
            oplot, [xright, xright], [0.0, 2.0], linestyle = 2, col = djs_icolor('green')
         wait, 1
         ;;IF islit EQ 25 THEN STOP 
      ENDIF
   ENDFOR
 ;;  old_slits = tset_slits
;;   width = old_slits[1].coeff[0, *]-old_slits[0].coeff[0, *]
;;   keep_slits = WHERE(width GE MINSLIT AND width LE MAXSLIT, nkeep)
;;   tset_proto = $
;;      { func    :    old_slits[0].FUNC, $
;;        xmin    :    old_slits[0].XMIN, $
;;        xmax    :    old_slits[0].XMAX, $
;;        coeff   :    dblarr(ncoeff, nkeep),  $
;;        dims    :    old_slits[0].DIMS, $
;;        xcorr_coeff : dblarr(ncoeff, nslit)  $
;;      }
;;   tset_slits = replicate(tset_proto, 2)
;;   tset_slits[0].COEFF[*, *] = old_slits[0].COEFF[*, keep_slits]
;;   tset_slits[1].COEFF[*, *] = old_slits[1].COEFF[*, keep_slits]
;;   tset_slits[0].XCORR_COEFF = old_slits[0].COEFF
;;   tset_slits[1].XCORR_COEFF = old_slits[1].COEFF
;; ;endelse
;;   nslit = nkeep
  
RETURN, tset_slits
END
