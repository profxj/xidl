FUNCTION TSPEC_WAVEIMG, piximg, ordermask, order_vec
  
  waveimg = double(0*piximg)
  ;; KLUDGE from P200 webpage
  ;wave_min = 1d4*[0.80778160D, 0.94296033D, 1.1297587D, 1.4103792D, 1.8760142D]
  ;wave_max = 1d4*[1.0629230D, 1.2379470D, 1.4835071D, 1.8501315D, 2.4644498D]
  wave_min = alog10(1d4*[0.80778160D, 0.94296033D, 1.1297587D, 1.4103792D, 1.8760142D])
  wave_max = alog10(1d4*[1.0629230D, 1.2379470D, 1.4835071D, 1.8501315D, 2.4644498D])
  dims = size(waveimg, /dim)
  ny = dims[1]
  
  norders = n_elements(order_vec)
  FOR iord = 0L, norders-1L DO BEGIN
     iorder = WHERE(ordermask EQ order_vec[iord])
     waveimg[iorder] = wave_min[iord] + double(piximg[iorder])*(wave_max[iord]-wave_min[iord])/double(ny-1)
  ENDFOR
  waveimg = 10.0d^waveimg
  RETURN, waveimg
END
