# REFL: Test the reflection algorithm
# Change CCD_GEOM to put angles in radians.  Review all real/double

include	<math.h>
include	"instrument.h"


define	D_1	20018.4D0	# (see PPLDIST)
# define	D_2	22215.5D0	# PPLDIST+86.5in=2197.1mm)
define	R_COLL	4394.2D0
define	K_COLL	-0.75
define	R_CURV	2124.71D0	# R_IMSURF here 83.65 in??
# define	D_1	20023.15D0
# define	D_2	22220.25D0
# define	R_CURV	2071.88D0

define	NCCD	8


procedure	t_refl()

double	tx, ty
double	t1, p1
double	t2, p2
double	t3, p3
real	xpix, ypix
double	x[3], r[3]		# variables for flens

real	clgetr()

begin
	tx = clgetr ("tx")
	ty = clgetr ("ty")
	t1 = clgetr ("t1")
	p1 = clgetr ("p1")
	t2 = clgetr ("t2")
	p2 = clgetr ("p2")
	t3 = clgetr ("t3")
	p3 = clgetr ("p3")
	xpix = clgetr ("xpix")
	ypix = clgetr ("ypix")

	call eprintf ("Inp angles: %7f %7f \n")
		call pargd (tx)
		call pargd (ty)

#	call setup0 (tx, ty, t1, p1, t2, p2, t3, p3, xpix, ypix)

#	call flens (tx, ty, x, r)	# tx, ty = location in mm in FP

end

procedure	setup0 (tx, ty, t1, p1, t2, p2, t3, p3, xpix, ypix)

double	tx, ty
double	t1, p1		# theta, phi of normal 1
double	t2, p2		# theta, phi of normal 2
double	t3, p3		# theta, phi of normal 3
real	xpix, ypix

real	ccd[NCCD,3]		# CCD geometrical params
double	r[3], a1[3,3], a2[3,3], a3[3,3], a4[3,3]

double	phi, theta
double	phin, thetan
double	tanx, tany
double	cost, sint, cosp, sinp

double	x, y
double	alpha, beta, gamma

real	ord, wmin, wmax
double	space, w
int	npts
int	i, n

int	ident_ccd()
int	clgeti()
real	clgetr()

begin
## TMP!
	n = 9
	call ccd_geom (ccd, 0., 0., 0.)

	ord = clgeti ("norder")
	space = 1000. / clgeti ("gmm")
	wmin = clgetr ("wmin")
	wmax = clgetr ("wmax")
	npts = clgeti ("npts")
	call printf ("# ord=%2f  space=%5f  tilt=%6.2f  tx=%6.3f  tx=%6.3f\n")
		call pargr (ord)
		call pargd (space)
		call pargd (t3)
		call pargd (tx)
		call pargd (ty)
# CALC x,y,z
	tx = DEGTORAD (tx)
	ty = DEGTORAD (ty)
	tanx = tan (tx)
	tany = tan (ty)
	theta = atan (sqrt (tanx*tanx + tany*tany))
	phi = atan2 (tany, tanx)

	cosp = cos (phi)
	sinp = sin (phi)
	cost = cos (theta)
	sint = sin (theta)

	r[1] = cosp * sint
	r[2] = sinp * sint
	r[3] = cost

call eprintf ("                       Into x,y,z: %7f %7f %7f\n")
call pargd (r[1])
call pargd (r[2])
call pargd (r[3])

# This is NOT exact -- tmp TMP
	x = r[1] * 20023.15
	y = r[2] * 20023.15
call eprintf ("x,y; %6f %6f\n")
call pargd (x)
call pargd (y)


# COLLIMATOR:
	call coll_angle (x, y, theta, phi)

	thetan = theta + DEGTORAD(t1)
	phin = phi + DEGTORAD(0.) + HALFPI
# The equations above are not fair -- there should be a transformation, based
# on (thetax,thetay) into new (theta,phi).  This will be true everywhere ...
#	phin = phi + DEGTORAD(p1) + HALFPI

	cosp = cos (phin)
	sinp = sin (phin)
	cost = cos (thetan)
	sint = sin (thetan)

#	 cosp		 sinp		0
#	-cost*sinp	 cost*cosp	sint
#	 sint*sinp	-sint*cosp	cost

	a1[1,1] = cosp
	a1[1,2] = sinp
	a1[1,3] = 0.
	a1[2,1] = -cost*sinp
	a1[2,2] = cost*cosp
	a1[2,3] = sint
	a1[3,1] = sint*sinp
	a1[3,2] = -sint*cosp
	a1[3,3] = cost

# TENT MIRROR:
# this mirror is OK to leave in del-theta,phi
	thetan = DEGTORAD(71.5 + t2)
	phin = DEGTORAD(90. + p2) + HALFPI
	cosp = cos (phin)
	sinp = sin (phin)
	cost = cos (thetan)
	sint = sin (thetan)

	a2[1,1] = cosp
	a2[1,2] = sinp
	a2[1,3] = 0.
	a2[2,1] = -cost*sinp
	a2[2,2] = cost*cosp
	a2[2,3] = sint
	a2[3,1] = sint*sinp
	a2[3,2] = -sint*cosp
	a2[3,3] = cost

# MIRROR:
# Better in thetax, thetay description (see above)
# for GRATING: need to add the third rotation
# Note the hack with phin + PI.  This is needed to keep the transformed
# x, y axes from flipping, wrt the original x,y.  Not a problem wrt reflections
# but if we want to work _within_ a system it is a problem. Cf. camera also
# note that this _forces_ us to work in a particular hemisphere, and thus we
# must make use of the negative theta as needed.

	thetan = DEGTORAD(t3)
	phin = DEGTORAD(p3) + HALFPI
	if (p3 > 0.) {
		phin = phin + PI
		thetan = -thetan
	}
	cosp = cos (phin)
	sinp = sin (phin)
	cost = cos (thetan)
	sint = sin (thetan)

	a3[1,1] = cosp
	a3[1,2] = sinp
	a3[1,3] = 0.
	a3[2,1] = -cost*sinp
	a3[2,2] = cost*cosp
	a3[2,3] = sint
	a3[3,1] = sint*sinp
	a3[3,2] = -sint*cosp
	a3[3,3] = cost

# CAMERA
# again, tranformation from theta-x,y is better, although this is not ridiculous
	thetan = DEGTORAD(2.33)
	phin = DEGTORAD(90.) + HALFPI
		phin = phin + PI
		thetan = -thetan
	cosp = cos (phin)
	sinp = sin (phin)
	cost = cos (thetan)
	sint = sin (thetan)

	a4[1,1] = cosp
	a4[1,2] = sinp
	a4[1,3] = 0.
	a4[2,1] = -cost*sinp
	a4[2,2] = cost*cosp
	a4[2,3] = sint
	a4[3,1] = sint*sinp
	a4[3,2] = -sint*cosp
	a4[3,3] = cost

# call eprintf ("A %7.4f %7.4f %7.4f\n  %7.4f %7.4f %7.4f\n  %7.4f %7.4f %7.4f\n")
# call pargd (a4[1,1])
# call pargd (a4[1,2])
# call pargd (a4[1,3])
# call pargd (a4[2,1])
# call pargd (a4[2,2])
# call pargd (a4[2,3])
# call pargd (a4[3,1])
# call pargd (a4[3,2])
# call pargd (a4[3,3])
	call refl (r, a1)
	call refl (r, a2)
#	call refl (r, a3)
# call eprintf ("Before the grating system    x,y,z: %7f %7f %7f\n")
# call pargd (r[1])
# call pargd (r[2])
# call pargd (r[3])
	call gen_xfm (r, a3, YES)

# call eprintf ("Within the grating system    x,y,z: %7f %7f %7f\n")
# call pargd (r[1])
# call pargd (r[2])
# call pargd (r[3])

	tx = atan2 (-r[1], -r[3])
	ty = atan2 (-r[2], -r[3])

## HERE IS THE GRATING EQUATION:
## m*wave = -n * space * cos(gamma) * ( sin(beta) + sin(alpha) )
	alpha = -ty
	gamma = -tx   # both alpha, gamma seem reversed; but gamma must be *-1
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))
call eprintf ("##* gamma = %7.3f\n")
call pargd (RADTODEG(gamma))

# print out normal (TMP)
	beta = 0.
	call grat_to_ics (beta, beta, 381.5D0, a4, a3, ccd, n, xpix, ypix)
	n = ident_ccd (xpix, ypix)
	call ics_to_ccd (xpix, ypix, ccd, n, xpix, ypix)

	call printf ("## %7.1f %7.1f  #  normal\n")
		call pargr (xpix)
		call pargr (ypix)

	do i = 0, npts-1 {
		w = wmin + (wmax-wmin)/(npts-1.)*i

		beta = asin ((ord * w / space / cos (gamma)) - sin (alpha))
		call grat_to_ics (beta, gamma, 381.5D0, a4, a3, xpix, ypix)
		n = ident_ccd (xpix, ypix)
		call ics_to_ccd (xpix, ypix, ccd, n, xpix, ypix)

		call printf ("%7.1f %7.1f  #  %7.1f \n")
			call pargr (xpix)
			call pargr (ypix)
			call pargd (w*10000.)
	}

end

procedure	refl (r, a)

double	r[3]
double	a[3,3]

double	rp[3]
int	i, j

begin
# transform
	do i = 1, 3 {
	    rp[i] = 0.
	    do j = 1, 3 {
		rp[i] = rp[i] + a[i,j] * r[j]
	    }
	}
		
# for a reflection, change (zp) --> (-zp)

	rp[3] = -rp[3]
	
# re-transform
	do i = 1, 3 {
	    r[i] = 0.
	    do j = 1, 3 {
		r[i] = r[i] + a[j,i] * rp[j]
	    }
	}
end

#
# GEN_XFM: general transform of r[3] into another CS desribed by a; "forward"
# is YES/NO to describe if xform is into or ou-of CS
# Note that the appropriate operation (eg transmission, reflection) must be
# applied afterward
#

procedure	gen_xfm (r, a, forward)

double	r[3]
double	a[3,3]
int	forward

double	rp[3]
int	i, j

begin
# transform
	if (forward == YES) {
	    do i = 1, 3 {
		rp[i] = 0.
		do j = 1, 3 {
			rp[i] = rp[i] + a[i,j] * r[j]
		}
	    }

	    do i = 1, 3
		r[i] = rp[i]

	} else {
		
	    do i = 1, 3
		rp[i] = r[i]

	    do i = 1, 3 {
		r[i] = 0.
		do j = 1, 3 {
			r[i] = r[i] + a[j,i] * rp[j]
		}
	    }
	}
end


procedure	coll_angle (xp, yp, sys, tc, pc)
# NOTE that returned tc, pc are the downward normal

double	xp, yp			# projected x,y in plane of telescope
double	sys[NPARAM]		# system parameters
double	tc, pc			# returned theta, phi of collimator surface

double	rp, rc			# radius projected, collimator
double	hm			# height at mask
double	cott			# cotangent theta
double	d2			# pupil to collimator dist
double	d, k			# work factors

begin
	pc = atan2 (yp, xp)		# phi will be same on coll

	rp = sqrt (xp*xp + yp*yp)
	hm = R_CURV - sqrt (R_CURV*R_CURV - rp*rp)

	d2 = D_1 + COL_DST(sys)

	cott = (D_1 - hm) / rp
	k = 1. + (1. + K_COLL) * cott*cott
	d = d2 * (1. + K_COLL)

# The following is general for conic sections.  In practice, a parabola is fine
# Note the switch in sign for quadratic root
	if (R_COLL - d > 0.) {
		rc = (R_COLL - d) / k *
	(sqrt (cott*cott + d2 * k * (2.*R_COLL - d) / (R_COLL - d)**2) - cott)
	} else {
		rc = (d - R_COLL) / k *
	(sqrt (cott*cott + d2 * k * (2.*R_COLL - d) / (R_COLL - d)**2) + cott)
	}
# This is for parabola:
#	rc = R_COLL * (sqrt (cott*cott + 2. * d2 / R_COLL ) - cott)
	
# call eprintf ("phi-c, rc: %5f %5f\n")
# call pargd (RADTODEG(pc))
# call pargd (rc)

# The general conic form (important)
	rc = rc / R_COLL
	tc = atan (rc / sqrt (1. - (1.+K_COLL)*rc*rc))
end


procedure	ics_to_grat (x, y, sys, a, ag, beta, gamma)

real	x, y
double	sys[NPARAM]			# system parameters
double	a[3,3]				# camera tranformation
double	ag[3,3]				# grating tranformation
double	beta, gamma			# returned angles


double	theta, phi
double	xp, yp, rp 
double	xpp, ypp
double	cosp, sinp, cost, sint

double	r[3]			# coord vector

begin

# 8. Remove the mosaic transform
	cosp = cos (-MOS_ROT(sys))
	sinp = sin (-MOS_ROT(sys))
	xpp = x - X_OPT(sys)
	ypp = y - Y_OPT(sys)
	xp =  cosp * xpp + sinp * ypp
	yp = -sinp * xpp + cosp * ypp

# 7. The display coords are different than the x,y,z convention adopted so far:

	xp = -xp
	yp = -yp

# 6. Convert to x,y in focal plane:
	rp = sqrt (xp*xp + yp*yp)
	theta = atan (rp * PIX_SZ / CAM_FOC(sys))
	phi = atan2 (yp, xp)

# 5. Camera distort theta (phi unchanged)
	call cam_distort (theta, theta, NO)

# 4. Convert into x,y,z:

	cosp = cos (phi)
	sinp = sin (phi)
	cost = cos (theta)
	sint = sin (theta)

	r[1] = cosp * sint
	r[2] = sinp * sint
	r[3] = cost

# 3. transform out of camera system
	call gen_xfm (r, a, NO)

# 2. transform into grating system
	call gen_xfm (r, ag, YES)

# 1. convert x,y,z into gamma, beta; note sign reversal of beta

	beta = -atan2 (r[2], r[3])		# NB the assoc. of y w/beta
	gamma = asin (r[1])			# NB the assoc. of x w/gamma!
	

call eprintf ("poss. sign err. in x,y\n")
call eprintf (" beta: %6f  gamma: %6f \n")
call pargd (RADTODEG(beta))
call pargd (RADTODEG(gamma))

end




procedure	grat_to_ics (beta, gamma, sys, a, ag, x, y)

double	beta, gamma			# returned angles
double	sys[NPARAM]			# system parameters
double	a[3,3]				# camera tranformation
double	ag[3,3]				# grating tranformation
real	x, y				# returned pixel values in ICS


double	theta, phi
double	xp, yp, rp
double	cosp, sinp
# double	cost, sint
# double	tanx, tany

double	r[3]			# coord vector

begin
# convert beta, gamma into x,y,z (cf Schroeder p259); note sign reversal of beta
	r[1] = sin (gamma)
	r[2] = sin (-beta) * cos (gamma)
	r[3] = cos (-beta) * cos (gamma)

# call eprintf ("x,y,z: %7f %7f %7f\n")
# call pargd (r[1]) ; call pargd (r[2]) ; call pargd (r[3])

# transform out of grating system
	call gen_xfm (r, ag, NO)

# call eprintf ("Out of grating system x,y,z: %7f %7f %7f\n")
# call pargd (r[1]) ; call pargd (r[2]) ; call pargd (r[3])

# transform into camera system
	call gen_xfm (r, a, YES)

# call eprintf ("In to camera system x,y,z: %7f %7f %7f\n")
# call pargd (r[1]) ; call pargd (r[2]) ; call pargd (r[3])

# Convert into theta, phi
	phi = atan2 (r[2], r[1])
	theta = atan2 (sqrt (r[1]*r[1] + r[2]*r[2]), r[3])

# Camera distort theta (phi unchanged)
	call cam_distort (theta, theta, YES)

# Convert to x,y in focal plane:

	rp = tan (theta) * CAM_FOC(sys) / PIX_SZ
	xp = rp * cos (phi)
	yp = rp * sin (phi)

# The display coords are different than the x,y,z convention adopted so far:

	xp = -xp
	yp = -yp

# Apply the mosaic transform
	cosp = cos (MOS_ROT(sys))
	sinp = sin (MOS_ROT(sys))
	x =  cosp * xp + sinp * yp + X_OPT(sys)
	y = -sinp * xp + cosp * yp + Y_OPT(sys)


# call eprintf ("    x,y in ICS: %5f %5f ")
# call pargr (x)
# call pargr (y)


end


procedure	ccd_to_ics (xccd, yccd, ccd, n, xics, yics)

real	xccd, yccd			# Input CCD values
real	ccd[NCCD,3]			# CCD geometrical params
int	n
real	xics, yics			# returned ICS values

double	cosa, sina, bx, by
double	x, y
double	xp, yp

begin
# TMP setup
	cosa = cos (ccd[n,3])
	sina = sin (ccd[n,3])
	bx = ccd[n,1]			# Refer to chip center
	by = ccd[n,2]			# Refer to chip center

# 8. Convert to Image CS:
	x = xccd - 1024.
	y = yccd - 2048.

	xp =  cosa * x - sina * y + bx
	yp =  sina * x + cosa * y + by

	xics = xp
	yics = yp

end

procedure	ics_to_ccd (xics, yics, ccd, n, xccd, yccd)

real	xics, yics			# Input ICS coords
real	ccd[NCCD,3]			# CCD geometrical params
int	n
real	xccd, yccd			# returned CCD values

double	cosa, sina, bx, by
double	xp, yp
double	x, y

begin
	cosa = cos (ccd[n,3])
	sina = sin (ccd[n,3])
	bx = ccd[n,1]			# These are center of chip
	by = ccd[n,2]			# These are center of chip

	xp = xics - bx
	yp = yics - by
	x =  cosa * xp + sina * yp + 1024.
	y = -sina * xp + cosa * yp + 2048.

# call eprintf ("  x,y in CCD system %d: %5f %5f \n")
# call pargi (n)
# call pargd (x)
# call pargd (y)

	xccd = x
	yccd = y

end

int	procedure ident_ccd (xim, yim)

real	xim, yim			# x,y in ICS (pix)

int	mx, my
real	xp, yp				# transformed
real	cosa, sina, xoff, yoff
real	ay, by
real	a1, b1, a2, b2, a3, b3

begin
# TMP: these will need to be input from somewhere ...
# NB REVIEW! Fix signs, etc
	cosa = cos (DEGTORAD(0.))
	sina = sin (DEGTORAD(0.))
	xoff = 0.
	yoff = 0.
	xp =  cosa * xim + sina * yim + xoff
	yp = -sina * xim + cosa * yim + yoff

# TMP: These boundaries currently faked
	ay = 0.
	by = 0.
	a1 = 0.
	a2 = 0.
	a3 = 0.
	b1 = -(2048. + 90.)		# was 66.7, now has dead Si space
	b2 = 0.
	b3 =  (2048. + 90.)
	
# what is y-height of dividing line?
	if (yp > (xp * ay + by))
		my = 2
	else
		my = 1

# find appropriate x-division
	if (xp < (a1 * yp + b1))
		mx = 1
	else if (xp < (a2 * yp + b2))
		mx = 2
	else if (xp < (a3 * yp + b3))
		mx = 3
	else
		mx = 4

	return ((my-1)*4+mx)
end
