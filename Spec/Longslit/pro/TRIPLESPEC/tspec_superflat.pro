function tspec_superflat, filenames, slitfile, pixfile $
                          , order_vec = order_vec, skyflat = SKYFLAT $
                          , NOSCAT = NOSCAT, ILLUM = ILLUMFINAL $
                          , TWEAK = TWEAK, TOL_EDG = TOL_EDG1 $
                          ,  EDG_TRIM = EDG_TRIM1 $
                          , CHK = CHK
    
  ;;scattfile = strtrim(archive_path, 2)+'scatt_template.fits'
  
  prog_name = 'tspec_superflat';

  nflat = n_elements(filenames)
  IF nflat LT 1 THEN message, 'No flat field files given in tspec_superflat'
  IF KEYWORD_SET(TOL_EDG1) THEN TOL_EDG = TOL_EDG1 $
  ELSE TOL_EDG = 3
  IF KEYWORD_SET(EDG_TRIM1) THEN EDG_TRIM=EDG_TRIM1 $
  ELSE EDG_TRIM=0

  ;; dist from edge for scattered light
  IF NOT KEYWORD_SET(scatt_shift) THEN scatt_shift = 2L 

  ;; create ordermask from tset_slits structure
  ;ordermask=xmrdfits(slitfile,0)
  tset_slits = xmrdfits(slitfile, 1)
  ordermask = tspec_ordermask(tset_slits) ;;, /fsr)
  dslit = tset_slits[1].COEFF[0, *]-tset_slits[0].COEFF[0, *]

  slitmask = long_slits2mask(tset_slits)
  slitmask_left  = long_slits2mask(tset_slits, xshift = -abs(scatt_shift))
  slitmask_right = long_slits2mask(tset_slits, xshift = abs(scatt_shift))
  scattmask = (slitmask_left EQ 0 AND slitmask_right EQ 0 AND slitmask EQ 0)

  ; This code snippet blocks out orders at the edge of the detector 
  ; that we don't use
  traceset2xy, tset_slits[0], yy, xx ; Left hand trace
  traceset2xy, tset_slits[1], yy2, xx2 ; Right hand trace
;  for yind=0, 2048-1 do begin
;     scattmask[0:(xx[yind,0]-25>1),yind]         = 0B
;     scattmask[(xx2[yind,20]+55<2047):2047,yind] = 0B
; endfor

  ;stop
  ximg = long_slits2x(tset_slits, edgmask = edgmask, TOL_EDG = TOL_EDG)
  piximg = xmrdfits(pixfile, 0)

  ;; these are the same for all images
  ;;scatt_model = xmrdfits(scattfile, 0)
    
  ;; which orders are we flattening
  IF NOT KEYWORD_SET(order_vec) THEN BEGIN 
     order_vec = reverse(dindgen(5) + 3)
     norders = 5L
  ENDIF ELSE norders = n_elements(order_vec)
  
  IF KEYWORD_SET(SKYFLAT) THEN specsamp = 1.5D $
  ELSE specsamp = 5L ;; spectral direction sampling in pixels
  spatsamp = 0.1D ;; spatial sampling for illum func in fractional order widths

  print, prog_name, ': Combining flat files using tspec_combine_flats.pro...'
  TSPEC_COMBINE_FLATS, FILENAMES=filenames, ORDERMASK=ordermask, edgmask=edgmask, $
                       IMAGE = image , OUTMASK = outmask ;, /ADD

  print, prog_name, ':    ...flat files combined.'
  dims = size(image, /dim)
  nx = dims[0]
  ny = dims[1]
  xarr = findgen(nx) # replicate(1.0, ny)

  IF NOT KEYWORD_SET(NOSCAT) THEN begin
     print, " "
     print, "Fitting 2D scattered light model for flat..."
     print, " "
;     stop
     scatt_model = tspec_scattlight(image, tset_slits, scatt_shift=scatt_shift)
  endif

  illum_flat = fltarr(nx, ny)
  scattpix = WHERE(scattmask AND image GT 1.0 AND image LT 5d4 AND $ ;; Bright flats exceed 5d4!!
                   (xarr GT 0 AND xarr LT (nx-1)))
  ;; fit for the amplitude of scattered light for this flat using
  ;; the spatial template
  IF NOT KEYWORD_SET(NOSCAT) THEN BEGIN
      scatt_amp = total(image[scattpix]*scatt_model[scattpix])/ $
        total(scatt_model[scattpix]^2)
      image = image - scatt_amp*scatt_model
   ENDIF
  
  ;; now normalize out the spectral and spatial depdendence
  FOR iorder = 0L, norders-1L DO BEGIN
      print, 'Generating flat field for order '+strtrim(order_vec[iorder], 2)
      inorder = where(ordermask EQ order_vec[iorder], nord)
      fitpix = where(ordermask  EQ order_vec[iorder]  $
                     AND finite(image)  $
                     AND abs(image) LE 6.0d4 $
                     AND NOT EDGMASK, nfit)
      if iorder EQ 4 then stop
      psort = sort(piximg[fitpix])
      med_spec_width = round(double(nfit)*double(specsamp)/double(ny))
      normspec1 = djs_median(image[fitpix[psort]], Width = med_spec_width $
                             , bound = 'reflect')
      ;; now smooth with a gausian 
      sig_res_spec = med_spec_width/15L
      nhalf =  long(sig_res_spec)*4L
      xkern = dindgen(2*nhalf+1)-nhalf
      kernel = gauss1(xkern, [0.0, sig_res_spec, 1.0])
      normspec = convol(normspec1, kernel, /edge_truncate)
      sset_spec = bspline_iterfit(piximg[fitpix[psort]], normspec $
                                  , upper = 5, lower = 5, yfit = yfit $
                                  , bkspace = 0.5d, maxiter = 20 $
                                  , maxrej = 10, /silent)
        
;      stop
      image[inorder] = image[inorder]/bspline_valu(piximg[inorder], sset_spec)
        
      ;; now fit out average illuimination function across the order
      fitpix2 = where(ordermask EQ order_vec[iorder]  $
                      AND finite(image)  $
                      AND abs(image) LE 6.0d4, nfit2)
        
      psort2 = sort(ximg[fitpix2])
      med_spat_width = round(double(nfit2)*double(spatsamp))
      normspat1 = djs_median(image[fitpix2[psort2]], width = med_spat_width $
                             , bound = 'reflect')
      ;; now smooth with a gausian 
      sig_res_spat = med_spat_width/15L
      nhalf =  long(sig_res_spat)*4L
      xkern = dindgen(2*nhalf+1)-nhalf
      kernel = gauss1(xkern, [0.0, sig_res_spat, 1.0])
      normspat = convol(normspat1, kernel, /edge_truncate)
      statinds = WHERE(ximg[fitpix2[psort2]] GT 0.1 AND $
                       ximg[fitpix2[psort2]] LT 0.9)
      moms = moment(normspat[statinds])
      mean = moms[0]
      normspat = normspat/mean
      nrm_flat_orig = image[fitpix2[psort2]]/mean
      sset_spat = bspline_iterfit(ximg[fitpix2[psort2]], normspat $
                                  , upper = 3, lower = 3, yfit = yfit2 $
                                  , bkspace = spatsamp/50.0d $
                                  , maxiter = 20, maxrej = 10 $
                                  , /silent)
      ;; don't allow super big corrections to protect edges
      normbsp = (bspline_valu(ximg[inorder], sset_spat)) < 2.0
      normbsp = normbsp > 0.5
;      stop
      image[inorder] = image[inorder]/normbsp
      illum_flat[inorder] = normbsp
      xleft_tweak = -1.0d
      xright_tweak = -1.0d
      IF (KEYWORD_SET(TWEAK)) THEN BEGIN
         ;; now tweak the slit position so that the slits cut off no less 
         ;; than 80% from the boundary
         step = 0.005
         inorm = WHERE(ximg[fitpix2[psort2]] GE 0.2 AND $
                       ximg[fitpix2[psort2]] LE 0.8)
         maxnorm = max(yfit2[inorm])
         FOR xleft = 0.5D, 0.0, -step DO BEGIN
            ynow = bspline_valu(xleft, sset_spat)
            IF ynow LE 0.9d*maxnorm THEN BEGIN 
               xleft_tweak = xleft
               BREAK
            ENDIF
         ENDFOR
         IF xleft_tweak NE -1.0 THEN BEGIN
            tset_slits[0].COEFF[0, iorder] = $
               tset_slits[0].COEFF[0, iorder] + xleft_tweak*dslit[iorder]
            splog, 'Tweaking position of LEFT boundary for order', 7L-iorder
         ENDIF
         FOR xright = 0.5D, 1.0, step DO BEGIN
            ynow = bspline_valu(xright, sset_spat)
            IF ynow LE 0.9d*maxnorm THEN BEGIN 
               xright_tweak = xright
               BREAK
            ENDIF
         ENDFOR
         IF xright_tweak NE -1.0 THEN BEGIN
            tset_slits[1].COEFF[0, iorder] = $
               tset_slits[1].COEFF[0, iorder] - $
               (1.0d - xright_tweak)*dslit[iorder]
           splog, 'Tweaking position of RIGHT boundary for order', 7L-iorder
         ENDIF  
      ENDIF
      IF KEYWORD_SET(CHK) THEN BEGIN
         yr = [0.2, 1.3]
         clr=getcolor(/load)
         msg = 'Illumination function fit, slit '+strtrim(iorder+1,2)+' of '+strtrim(norders,2)
         plot, ximg[fitpix2[psort2]], nrm_flat_orig $
               , psym = 3, yr = yr, xr = [-0.1, 1.1], /xsty, /ysty $
               , charthick = 2.0, thick = 3.0, title=msg
         oplot, ximg[fitpix2[psort2]], yfit2, col = clr.red $
                ,  thick = 2.0
         IF xleft_tweak NE -1.0 THEN $
            oplot, [xleft_tweak, xleft_tweak], yr, linestyle = 2 $
                   , col = clr.green $
                   , thick = 2.0
         IF xright_tweak NE -1.0 THEN $
            oplot, [xright_tweak, xright_tweak], yr, linestyle = 2 $
                   , col = clr.green $
                   , thick = 2.0
      ENDIF
      
   ENDFOR                       ; loop through orders

  illumfinal = illum_flat

  ;; Write tweaked slit boundary to file, overwriting previous
  IF KEYWORD_SET(TWEAK) THEN BEGIN
     ordermask = tspec_ordermask(tset_slits) ;;, /fsr) 
      mwrfits, ordermask, slitfile, /create
      mwrfits, tset_slits, slitfile
      RETURN,0  ;; Return tweaked slit boundaries
  ENDIF
  IF KEYWORD_SET(CHK) THEN BEGIN 
      thismask = fltarr(nx, ny)
      FOR kk = 0L, norders-1L DO BEGIN
         ii = where(ordermask EQ order_vec[kk] AND outmask)
         thismask[ii] = 1L
      ENDFOR
;      xatv, thismask*image, min = 0.9, max = 1.1, /block
  ENDIF

  return, image
  
  print, prog_name, ':  Complete.'

END
  


