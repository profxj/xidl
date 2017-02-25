FUNCTION GMOS_SLITS2MASK, tset_slits1, tset_slits2, tset_slits3, ystruct
  
  nx = tset_slits1[0].dims[0]
  ny = tset_slits1[0].dims[1]
  idim = size(tset_slits1[0].COEFF, /dim)
  nslit = idim[1]
  gmos_gaps, nx, ny, ygap1 = ygap1, ygap2 = ygap2
  nmy = 3*ny + ygap1 + ygap2
  
  gmos_traceset3xy, tset_slits1, tset_slits2, tset_slits3 $
                    , left = xx1, right = xx2, yrow = yrow, NSLIT = NSLIT
  
  ;; Generate the mask image
  slitmask = lonarr(nx, nmy)
  slitmask_temp = lonarr(nx, nmy)
  for islit = 0L, nslit-1L do begin
     for iy = 0L, nmy-1L do begin
        x1 = round(xx1[iy, islit])
        x2 = round(xx2[iy, islit])
         if (x1 GE x2) AND KEYWORD_SET(VERBOSE) THEN $
            splog, 'WARNING: Slit start and end positions appear to cross!'
         x1 = x1 > 0
         x2 = x2 < (nx-1)
         if (x1 LE x2) then begin
            if (total(slitmask_temp[x1:x2, iy]) GT 0) $
               AND KEYWORD_SET(VERBOSE) then $
                  splog, 'WARNING: Slits appear to overlap!'
            slitmask_temp[x1:x2, iy] = islit+1
         endif
      endfor
  endfor
  
  yarr = findgen(nmy) ## replicate(1.0, nx)
  FOR islit = 0L, nslit-1L DO BEGIN
     ipix = WHERE(slitmask_temp EQ (islit+1L) AND $
                  yarr GE ystruct[islit].YSLIT_MIN AND $
                  yarr LE ystruct[islit].YSLIT_MAX)
     slitmask[ipix] = slitmask_temp[ipix]
  ENDFOR
  RETURN, slitmask
END
