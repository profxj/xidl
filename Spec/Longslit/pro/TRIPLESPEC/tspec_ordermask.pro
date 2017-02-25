FUNCTION TSPEC_ORDERMASK, tset_slits, order_vec = order_vec
  nx = tset_slits[0].dims[0]
  ny = tset_slits[0].dims[1]
  yarr = replicate(1.0, nx) # findgen(ny)
  ordermask = long_slits2mask(tset_slits)
  ordermask[WHERE(ordermask GT 0)] = -ordermask[WHERE(ordermask GT 0)] + 8L
  order_vec = [7L, 6L, 5L, 4L, 3L]
  i7 = WHERE(ordermask EQ 7 AND yarr LT 1100)
  ordermask[i7] = 0L
  RETURN, ordermask
END
