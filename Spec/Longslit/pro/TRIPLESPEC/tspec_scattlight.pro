FUNCTION tspec_scattlight, image1, tset_slits, scatt_shift = scatt_shift $
                           , sig_res = sig_res $
                           , CHK = CHK, SCATT_NORM = SCATT_NORM
  
  ;slitfile = '/Users/jhennawi/mage/mage_archive_orders.fits'
  ;tset_slits = mrdfits(slitfile, 1)
  ;mage_proc, '/Users/jhennawi/mage/nite1/mage0035.fits', image1, hdr = hdr
  ;slitfile = '/Users/jhennawi/mage/mage_archive_orders.fits'
  ;tset_slits = mrdfits(slitfile, 1)
  ;tset_slits = mage_traceorders(image1)

  IF NOT KEYWORD_SET(scatt_shift) THEN scatt_shift = 10L
  dims = size(image1, /dim)
  nx = dims[0]
  ny = dims[1]
  scatt = fltarr(nx, ny)
  ;; allocate some things needed for scattered light removal
  xarr = findgen(nx) # replicate(1.0, ny)
  yarr = findgen(ny)## replicate(1.0, nx)
  IF NOT KEYWORD_SET(SIG_RES) THEN sig_res = 5L
  nres = 4*sig_res +1L
  x = (-nres/2 + findgen(nres)) # replicate(1.0, nres)
  y = (-nres/2 + findgen(nres)) ## replicate(1.0, nres)
  ygauss = gauss2(x, y, [0.0, 0.0, sig_res, sig_res, 1.0])
  kernel = ygauss/total(ygauss)
  
  slitmask_left  = long_slits2mask(tset_slits, xshift = -abs(scatt_shift))
  slitmask_right = long_slits2mask(tset_slits, xshift = abs(scatt_shift))
  slitmask_middle = long_slits2mask(tset_slits)

  scattmask = (slitmask_left EQ 0 AND slitmask_right EQ 0 AND slitmask_middle EQ 0)
  traceset2xy, tset_slits[0], yy, xx ; Left hand trace
  traceset2xy, tset_slits[1], yy2, xx2 ; Right hand trace
  ;for yind = 0, 2048-1 do begin
  ;   scattmask[0:(xx[yind,0]-25>0),yind]     = 0B
  ;   scattmask[(xx2[yind,20]+55<2047):2047,yind] = 0B
  ;endfor

  ;; interpolate over any bad pixels so as not to break fits
  image2 = djs_maskinterp(image1, ((image1 LT 1.0 OR image1 GT 2d4) $
                          AND scattmask), iaxis = 0)
  image  = djs_maskinterp(image2, ((image2 LT 1.0 OR image2 GT 2d4) $
                          AND scattmask), iaxis = 0)
  mask = image GT 1.0 AND image LT 2d4 AND (xarr GT 0 AND xarr LT (nx-1))
  
  med_val = median(image[where(scattmask AND mask)])
  xy2traceset, xarr, image, tset $
               , invvar = double(scattmask AND mask)/3.0/med_val $
               , func = 'chebyshev' $
               , xmin = 0.0, xmax = float(nx-1L) $
               , ncoeff = 5, yfit = slitfit, upper = 3, lower = 3, /silent
  xy2traceset, transpose(yarr), transpose(slitfit), yset $
               , func = 'chebyshev' $
               , xmin = 0.0, xmax = float(ny-1L) $
               , ncoeff = 5, yfit = yslitfit, upper = 3, lower = 3, /silent
  smooth_set = yset
  nhalf =  20L
  xkern = dindgen(2*nhalf+1)-nhalf
  kern1 = gauss1(xkern, [0.0, 4.0, 1.0])
  FOR k = 0, n_elements(tset.COEFF[*, 0])-1L DO $
    smooth_set.COEFF[k, *] = convol(transpose(yset.coeff[k, *]), kern1 $
                                    , /edge_truncate)
;djs_median(transpose(tset.coeff[k, *]) $
                              ;           ,  width = 7, bound = 'reflect')
  traceset2xy, smooth_set, yrows, scattsmooth
  img_scatt = convol(transpose(scattsmooth), kernel, /edge_truncate)
;  img_scatt = transpose(scattsmooth)
  mom = moment(img_scatt)
  scatt_norm = mom[0]
  RETURN, img_scatt
END
