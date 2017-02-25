FUNCTION LUCI_EXTRACT_ARC, arcimg, arcivar, trace, box_rad = box_rad $
                           , var = sub_var

  if NOT keyword_set(box_rad) then box_rad = 5L

  dims = size(arcimg, /dim)
  nx = dims[0]
  ny = dims[1]
  arc1d = fltarr(ny)
  FOR j = 0L, ny-1L DO BEGIN
     left  = floor(trace[j] - BOX_RAD)
     right = ceil(trace[j] + BOX_RAD)
     sub_arc  = arcimg[left:right, j]
     sub_ivar = arcivar[left:right, j] > 0.
     sub_var = 1.0/(sub_ivar + (sub_ivar EQ 0.0))
     djs_iterstat, sub_arc, invvar = sub_ivar $
                   , mean = mean, sigrej = 3.0, mask = mask, median = median
     arc1d[j] = mean
     IF total(mask) NE 0 THEN var = total(sub_var*mask)/total(mask)^2 $
     ELSE var = total(sub_var)
  ENDFOR

  RETURN, arc1d
END
