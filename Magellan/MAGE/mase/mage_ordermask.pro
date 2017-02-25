FUNCTION MAGE_ORDERMASK, tset_slits
  nx = tset_slits[0].dims[0]
  ny = tset_slits[0].dims[1]
  yarr = replicate(1.0, nx) # findgen(ny)
  ordermask = long_slits2mask(tset_slits)
  ordermask[WHERE(ordermask GT 0)] = -ordermask[WHERE(ordermask GT 0)] + 21L
  ;; cut off order 19,20 where they end
  ord20mask = (ordermask NE 20) OR (yarr GT 850.0 AND ordermask EQ 20)
  ord19mask = (ordermask NE 19) OR (yarr GT 275.0 AND ordermask EQ 19)
  ;; cut off order 6 for lam > 10300A 
  ord6mask  = (ordermask NE 6)  OR (yarr LT 1020 AND ordermask EQ 6)
  ordermask = ordermask*ord20mask*ord19mask*ord6mask
  RETURN, ordermask
END
