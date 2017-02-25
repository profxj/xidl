PRO GMOS_TRACESET3XY, tset_slits1, tset_slits2, tset_slits3 $
                      , LEFT = LEFT, RIGHT = RIGHT, YROW = YROW $
                      , NSLIT = NSLIT

  nx = tset_slits1[0].DIMS[0]
  ny = tset_slits1[0].DIMS[1]
  idim = size(tset_slits1[0].COEFF, /dim)
  nslit = idim[1]

  gmos_gaps, nx, ny, ygap1 = ygap1, ygap2 = ygap2
  nmy = 3*ny + ygap1 + ygap2

  left = fltarr(nmy, nslit)
  right = fltarr(nmy, nslit)
  
  traceset2xy, tset_slits1[0], rows, left1
  traceset2xy, tset_slits1[1], rows, right1
  
  traceset2xy, tset_slits2[0], rows, left2
  traceset2xy, tset_slits2[1], rows, right2
  
  traceset2xy, tset_slits3[0], rows, left3
  traceset2xy, tset_slits3[1], rows, right3
  
  left[0:ny-1L, *] = left3
  left[ny + ygap1:2*ny + ygap1-1L, *] = left2
  left[2*ny+ygap1+ygap2:*, *] = left1
 
  right[0:ny-1L, *] = right3
  right[ny + ygap1:2*ny + ygap1-1L, *] = right2
  right[2*ny+ygap1+ygap2:*, *] = right1
 
  ;; Now interpolate in the gaps: ygap1
  xvec = [lindgen(ny), ny + ygap1 + lindgen(ny)]
  xint1 = ny + lindgen(ygap1)
  FOR islit = 0L, nslit-1L DO BEGIN
     left[ny:ny+ygap1-1L, islit] = $
        interpol([left3[*, islit], left2[*, islit]], xvec, xint1)
     right[ny:ny+ygap1-1L, islit] = $
        interpol([right3[*, islit], right2[*, islit]], xvec, xint1)
  ENDFOR
  ;; Now interpolate in the gaps: ygap2
  xvec = [ny + ygap1 + lindgen(ny), 2*ny + ygap1 + ygap2 + lindgen(ny)]
  xint2 = 2*ny + ygap1 + lindgen(ygap2)
  FOR islit = 0L, nslit-1L DO BEGIN
     left[2*ny+ygap1:2*ny+ygap1+ygap2-1L, islit] = $
        interpol([left2[*, islit], left1[*, islit]], xvec, xint2)
     right[2*ny+ygap1:2*ny+ygap1+ygap2-1L, islit] = $
        interpol([right2[*, islit], right1[*, islit]], xvec, xint2)
  ENDFOR
  YROW = findgen(nmy) # replicate(1.0, nslit)
  
  RETURN
END
