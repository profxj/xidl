pro mage_makeflat, path=flatpath, blue=bluefiles, red=redfiles $
                   , green = greenfiles, orders = orderfile $
                   , arcfile = arcfile, piximg = piximg, illum = illumfiles $
                   , output = outfile, CHK = CHK, ALLCHK = ALLCHK $
                   , SKYFLATRED = SKYFLATRED, NO_ILLUM=no_illum
  norders=15
  
  ;ported from mage_makepix... replaces the old Pixel image tab
  tset_slits = xmrdfits(orderfile, 1)
  mage_proc, arcfile, arcimg
  
  pixset = long_wavepix(arcimg, tset_slits, fwhm=3.0 $
                        , box_radius = 3.0 $
                        , sig_thresh = 3.0 $
                        , pkwdth = 5.0 $
                        , TOLER = 2.0 $
                        , CHK = ARCCHK, /VERBOSE)

  piximg = long_wpix2image(pixset, tset_slits)
  pixfile = 'piximg_flats.fits'
  mwrfits, piximg, pixfile, /create
  mwrfits, pixset, pixfile
   
  
  
  IF NOT FILE_TEST('Flat', /DIR) THEN FILE_MKDIR, 'Flat'

  ;; Run with an illumflat and tweak slit boundaries
  if not keyword_set(NO_ILLUM) then begin
     input = illumfiles
     junkflat = mage_superflat(input[0], orderfile, pixfile $
                               , /skyflat, /tweak, chk=allchk)
  endif
  
  ;;Create slit illumination function
  if not keyword_set(NO_ILLUM) then begin
     input = illumfiles
     junkflat = mage_superflat(input, orderfile, pixfile $
                               , /skyflat, illum = illumflat, chk=allchk)
  endif else begin
     illumflat = piximg
     illumflat[*] = 1.
  endelse
     
  input  =  bluefiles
  blueflat   = mage_superflat(input, orderfile, pixfile $
                              , order_vec = [20, 19, 18] $
                              , /skyflat, chk=allchk)
  input =  greenfiles
  greenflat  = mage_superflat(input, orderfile, pixfile $
                              , order_vec = [17, 16, 15, 14] $
                              , /skyflat, chk=allchk)
  input   =  redfiles
  redflat    = mage_superflat(input, orderfile, pixfile $
                              , order_vec = [13, 12, 11, 10, 9, 8, 7, 6] $
                              , skyflat = SKYFLATRED, chk = allchk)

  EDG_TRIM=[4L,4L]
  ;; Remove EDG_TRIM pixels from each side. 
  IF KEYWORD_SET(EDG_TRIM) THEN BEGIN
     tset_slits=mrdfits(orderfile,1)
     FOR iorder = 0L, norders-1L DO BEGIN
        tset_slits[0].COEFF[0, iorder] = $
           tset_slits[0].COEFF[0, iorder] + EDG_TRIM[0]
        tset_slits[1].COEFF[0, iorder] = $
           tset_slits[1].COEFF[0, iorder] - EDG_TRIM[1]
     ENDFOR
     ordermask=mage_ordermask(tset_slits) 
     mwrfits, ordermask, orderfile, /create
     mwrfits, tset_slits, orderfile
  ENDIF
  tset_slits=mrdfits(orderfile,1)
  ordermask=mage_ordermask(tset_slits) 
  
  ;;TOL_EDG=[0]
  ;;ximg = long_slits2x(tset_slits, edgmask = edgmask, TOL_EDG = TOL_EDG $
  ;;                    , nslit = norders)
  ordermask=mage_ordermask(tset_slits)
  flat = blueflat * (ordermask LE 19 AND ordermask GE 18) + $
         greenflat * (ordermask LE 17 AND ordermask GE 14) + $
         redflat * (ordermask LE 13 AND ordermask GE 6)
  ;; force order 20 to be 1 everywhere  
  flat[where(ordermask EQ 20)] = 1.0
  illumflat[where(ordermask EQ 20)]=1.0
  
  ;;unit1=(edgmask OR flat GT 3.0 OR flat LE 0.0 OR ordermask EQ 0)
  unit1=(flat GT 3.0 OR flat LE 0.0 OR ordermask EQ 0)
  unitind1 = WHERE(unit1, nunit1)
  IF nunit1 GT 0 THEN flat[unitind1] = 1.0d
  
  ;;unit2=(edgmask OR illumflat GT 3.0 OR illumflat LE 0.2 OR ordermask EQ 0)
  unit2=(illumflat GT 3.0 OR illumflat LE 0.2 OR ordermask EQ 0)
  unitind2 = WHERE(unit2, nunit2)
  IF nunit2 GT 0 THEN illumflat[unitind2] = 1.0d

  IF KEYWORD_SET(CHK) THEN BEGIN
     xatv, flat, min = 0.9, max = 1.1, /block
     xatv, illumflat, /block, min = 0.5, max = 1.5
  ENDIF
  hdr=xheadfits(illumfiles[0])
  slit=sxpar(hdr,'SLITNAME')
  mwrfits, flat, 'Flat/Pixflat.fits', /create
  mwrfits, illumflat, 'Flat/Illumflat_' + slit + '.fits', /create
  
end
