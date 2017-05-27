pro ics_to_ccd,xics, yics, ccd, n, xccd, yccd

; see refl.x

; Converts ICS coordinates to CCD system

cosa=cos(ccd[n,2])
sina=sin(ccd[n,2])

bx=ccd[n,0]   ; center of chip coords
by=ccd[n,1]

	xp = xics - bx
	yp = yics - by
	x =  cosa * xp + sina * yp + 1024.
	y = -sina * xp + cosa * yp + 2048.

	xccd = x
	yccd = y
return
end
