function xgety_plt, ytmp, pos, xymnx, size, STRCT=strct

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'y = xgety_plt(ypix,pos,xymnx,size)'
    print,'or y = xgety_plt(strct, /STRCT)'
    return, -1
  endif 

; Allow for structure

  if keyword_set( STRCT ) then begin
      ypix = ytmp.ycurs
      pos = ytmp.pos
      xymnx = ytmp.xymnx
      idx = tag_exist(ytmp, 'winsize', /top_level)
      if idx EQ 0 then size = ytmp.size else size = ytmp.winsize
  endif else ypix = ytmp
  
; Returns the data value from the pixel value

  m = (xymnx[3] - xymnx[1]) / (size[1] * (pos[3] - pos[1]))
  b = xymnx[3] - m * size[1] * pos[3]
  return, m*ypix + b
end

