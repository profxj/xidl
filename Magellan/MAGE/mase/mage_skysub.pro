function mage_skymodel, sciimg, sciivar, piximg=piximg, orderimg=orderimg, $
                 tset_slits=tset_slits

  bsp  = 0.6
  FWHM = 3.0

  objstruct1 = long_objfind(sciimg, tset_slits = tset_slits, invvar = sciivar $
                          , skymask = skymask, objmask = objmask $
                          , nperslit = 1L, peakthresh = reduxthresh $
                          , fwhm = FWHM, ISLIT = ISLIT)

  slitmask = long_slits2mask(tset_slits)
  ximg = long_slits2x(tset_slits, edgmask = edgmask)

  skyimage = long_skysub(sciimg, sciivar, piximg, slitmask, skymask, edgmask $
                       , bsp = bsp, ISLIT = ISLIT)

  return, skyimage

end
