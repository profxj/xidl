function xgetx_plt, xtmp, pos, xymnx, size, STRCT=strct

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x = xgetx_plt(xpix,pos,xymnx,size)'
    print,'or Syntax - ' + $
             'x = xgetx_plt(strct, /STRCT)'
    return, -1
  endif 

; Allow for structure

  if keyword_set( STRCT ) then begin
      xpix = xtmp.xcurs
      pos = xtmp.pos
      xymnx = xtmp.xymnx
      idx = tag_exist(xtmp, 'winsize', /top_level)
      if idx EQ 0 then size = xtmp.size else size = xtmp.winsize
  endif else xpix = xtmp
  
; Returns the data value from the pixel value

  m = (xymnx[2] - xymnx[0]) / (size[0] * (pos[2] - pos[0]))
  b = xymnx[2] - m * size[0] * pos[2]
  return, m*xpix + b
end
