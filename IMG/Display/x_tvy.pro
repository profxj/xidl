function x_tvy, tv, INTG=intg

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'y = x_tvy(tv, /INTG) '
    return, -1
  endif 

; Returns the data value from the pixel value

  m = (tv.xymnx[3] - tv.xymnx[1]) / (tv.winsize[1] * (tv.pos[3] - tv.pos[1]))
  b = tv.xymnx[3] - m * tv.winsize[1] * tv.pos[3]

  if keyword_set( INTG ) then return, nint(m*tv.ycurs + b) $
    else return, m*tv.ycurs + b 

end
