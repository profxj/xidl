pro match_old, trans, x, y, oldx, oldy, offx, offy

	trans = intarr(n_elements(x))
	tmpx = oldx + offx
	tmpy = oldy + offy

	for i=0,n_elements(x)-1 do begin
	  max = 9999.
	  for j=0,n_elements(oldx)-1 do begin
	    dist = sqrt( (x(i)-tmpx(j))^2 + (y(i)-tmpy(j))^2 )
	    if(dist LT max) then begin
		max = dist
	 	svj = j
	    endif
          endfor
	  if(max LT 1.5) then trans(i) = svj $
	  else trans(i) = 9999
	endfor

return
end
	  
