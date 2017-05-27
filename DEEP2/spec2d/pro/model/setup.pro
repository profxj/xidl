pro setup,e1,a2,a3,a4,ccd,sys

; see trace.x

; COLL ERROR (first order):
; In this case, we must remove the phi we put in, hence the more complex form

HALFPI=!dpi/2.



thetan = 2.*sys.COL_ERR
phin   = sys.COL_PHI + HALFPI

cosp   = cos(phin)
sinp   = sin(phin)
cost   = cos(thetan)
sint   = sin(thetan)

e1=dblarr(3,3)

e1[0,0] =1-(1-cost)*sinp*sinp
e1[0,1] =  (1-cost) * cosp*sinp
e1[0,2] = -sinp*sint
e1[1,0] =  (1-cost) * cosp*sinp
e1[1,1] =  1 - (1-cost)*cosp*cosp	# sinp*sinp + cost*cosp*cosp
e1[1,2] =  cosp*sint
e1[2,0] =  sint*sinp
e1[2,1] = -sint*cosp
e1[2,2] =  cost



; TENT MIRROR:
; this mirror is OK to leave in del-theta,phi

thetan = sys.TNT_ANG
phin = sys.TNT_PHI + HALFPI

cosp   = cos(phin)
sinp   = sin(phin)
cost   = cos(thetan)
sint   = sin(thetan)

a2=dblarr(3,3)
a2[0,0] = cosp
a2[0,1] = sinp
a2[0,2] = 0.
a2[1,0] = -cost*sinp
a2[1,1] = cost*cosp
a2[1,2] = sint
a2[2,0] = sint*sinp
a2[2,1] = -sint*cosp
a2[2,2] = cost


; MIRROR/GRATING:
; Better in thetax, thetay description (see above)
; for GRATING: need to add the third rotation
; Note the hack with phin + PI.  This is needed to keep the transformed
; x, y axes from flipping, wrt the original x,y.  Not a problem wrt reflections
; but if we want to work _within_ a system it is a problem. Cf. camera also
; note that this _forces_ us to work in a particular hemisphere, and thus we
; must make use of the negative theta as needed.

thetan 	= -sys.MU
xsi    	= sys.GR_ZERR
rhon	= sys.GR_YERR

cost = cos (thetan)
sint = sin (thetan)
cosx = cos (xsi)
sinx = sin (xsi)
cosr = cos (rhon)
sinr = sin (rhon)

; ... below assumes phin=0. (ie adopts the roll/yaw approach)
a3=dblarr(3,3)

a3[0,0] =  cosx*cosr
a3[0,1] =  sint*sinr + cost*sinx*cosr
a3[0,2] = -cost*sinr + sint*sinx*cosr
a3[1,0] = -sinx	
a3[1,1] =  cost*cosx
a3[1,2] =  sint*cosx
a3[2,0] =  cosx*sinr
a3[2,1] = -sint*cosr + cost*sinx*sinr
a3[2,2] =  cost*cosr + sint*sinx*sinr


; CAMERA
; again, tranformation from theta-x,y is better, although this is not ridiculous
	
thetan 	= sys.CAM_ANG
phin 	= sys.CAM_PHI + HALFPI
phin 	= phin + !dpi
thetan 	= -thetan
cosp = cos (phin)
sinp = sin (phin)
cost = cos (thetan)
sint = sin (thetan)

a4=dblarr(3,3)
a4[0,0] = cosp
a4[0,1] = sinp
a4[0,2] = 0.
a4[1,0] = -cost*sinp
a4[1,1] = cost*cosp
a4[1,2] = sint
a4[2,0] = sint*sinp
a4[2,1] = -sint*cosp
a4[2,2] = cost

; CCD: define the geometry of the mosaic
	
ccd=ccd_geom(sys)

return
end
