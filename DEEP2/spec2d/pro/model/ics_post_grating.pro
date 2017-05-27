pro ics_post_grating,x,y,a,sys,r

; see qtrmap.x


;defsysv,'!PIX_SZ',0.015d0

npts=n_elements(x)
x=reform(double(x),npts)
y=reform(double(y),npts)


; a = camera transformation
; r = coord vector

; 8. Remove the mosaic transform
	cosp = cos (-sys.MOS_ROT)
	sinp = sin (-sys.MOS_ROT)
	xpp = x - sys.X_OPT
	ypp = y - sys.Y_OPT
	xp =  cosp * xpp + sinp * ypp
	yp = -sinp * xpp + cosp * ypp


; 7.The display coords are different than the x,y,z convention adopted so far:

	xp = -xp
	yp = -yp

; 6. Convert to x,y in focal plane:
	rp = sqrt (xp*xp + yp*yp)
	theta = atan(rp * !PIX_SZ / sys.CAM_FOC)
	phi = atan(yp, xp)

; 5. Camera distort theta (phi unchanged)
	
	cam_distort,theta, theta

; 4. Convert into x,y,z:

	cosp = cos(phi)
	sinp = sin(phi)
	cost = cos(theta)
	sint = sin (theta)

	r=dblarr(npts,3)
	r[*,0]=cosp*sint
	r[*,1]=sinp*sint
	r[*,2]=cost

;3.  transform out of camera system
	for i=0,npts-1 do begin
		rtmp=reform(r[i,*],3)
		gen_xfm,rtmp,a
		r[i,*]=rtmp
	endfor
return
end









