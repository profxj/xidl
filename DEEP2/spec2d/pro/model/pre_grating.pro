pro pre_grating,xmm,ymm,e1,a2,sys,r

; see qextr.x

;define	D_1	20018.4D0	; see deimos.h PPLDIST
;define	R_CURV	2124.71D0	; ditto; r=83.65in in this case = R_IMSURF
defsysv,'!D_1',20018.4D0
defsysv,'!R_CURV',2124.71D0

a1=dblarr(3,3)     
 
npts=n_elements(xmm)
pad=0.d0
xd=reform(double(xmm),npts)
yd=reform(double(ymm),npts)

; mask_to_proj already vectorized!
mask_to_proj,xd,yd,pad,x,y,pa

; convert to r[3]
rp=sqrt(x*x+y*y)
hm=!R_CURV-sqrt(!R_CURV*!R_CURV-rp*rp)
theta=atan (rp / (!D_1-hm))
phi = atan(y, x)

cosp = cos (phi)
sinp = sin (phi)
cost = cos (theta)
sint = sin (theta)

r=dblarr(npts,3)

r[*,0]=cosp*sint
r[*,1]=sinp*sint
r[*,2]=cost


; COLLIMATOR

; coll_angle already vectorized!
coll_angle,x,y,sys,theta,phi

for i=0,npts-1 do begin
	thetan = theta[i]
	phin   = phi[i]+!dpi/2.

	cosp = cos (phin)
	sinp = sin (phin)
	cost = cos (thetan)
	sint = sin (thetan)

	a1[0,0] = cosp
	a1[0,1] = sinp
	a1[0,2] = 0.
	a1[1,0] = -cost*sinp
	a1[1,1] = cost*cosp
	a1[1,2] = sint
	a1[2,0] = sint*sinp
	a1[2,1] = -sint*cosp
	a1[2,2] = cost
	rtmp=reform(r[i,*],3)
	refl,rtmp,a1
	gen_xfm,rtmp,e1,/forward
	refl,rtmp,a2
	r[i,*]=rtmp
endfor

return
end
