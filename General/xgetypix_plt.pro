function xgetypix_plt, y, pos, xymnx, size

  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'ypix = xgetypix_plt(y,pos,xymnx,size)'
    return, -1
  endif 

; Returns the data value from the pixel value

  m = size[1] * (pos[3] - pos[1]) / (xymnx[3] - xymnx[1]) 
  b = size[1] * pos[3] - m * xymnx[3]
  return, m*y + b
end
