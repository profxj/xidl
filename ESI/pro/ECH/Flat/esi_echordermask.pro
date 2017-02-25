FUNCTION ESI_ECHORDERMASK, tset_slits
  nx = tset_slits[0].dims[0]
  ny = tset_slits[0].dims[1]
  yarr = replicate(1.0, nx) # findgen(ny)
  ordermask = long_slits2mask(tset_slits)
  ordermask[WHERE(ordermask GT 0)] = -ordermask[WHERE(ordermask GT 0)] + 16L
  ;; cut off order 6 where it ends
  ord6mask = (ordermask NE 6) OR (yarr LT 2000.0 AND ordermask EQ 6)
  ord15mask = (ordermask NE 15) OR (yarr LT 3800 AND ordermask EQ 15)
  ordermask = ordermask*ord6mask*ord15mask
  RETURN, ordermask
END
