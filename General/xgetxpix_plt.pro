; Returns the pixel value from the data value
function xgetxpix_plt, x, pos, xymnx, size

  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'xpix = xgetxpix_plt(x,pos,xymnx,size)'
    return, -1
  endif 


  m = size[0] * (pos[2] - pos[0]) / (xymnx[2] - xymnx[0]) 
  b = size[0] * pos[2] - m * xymnx[2]
  return, m*x + b
end

