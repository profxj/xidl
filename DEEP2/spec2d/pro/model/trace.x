define	PRECOL	0		# Can be removed now ...
# define	PRECOL	48
# define	PRECOL	96

# TRACE: Trace mask points through DEIMOS
# NOTES: redefined p3 & GRATING XFM currently redefined for pure roll, but
# eventually probably want to return to old system as phi3 is probably a
# constant to be solved for

# TBD: there are 3 definitions for grating xfm in setup -- trim [commented out]
# TBD: implement errors in appropriate places
# TBD: resolve CCD description (done?)
# TBD: resolve p3/roll3 grating issue (see below)
# TBD: resolve error in xfmx: DONE, but more compact way to pass xfms desired ..
# TBD: put in "off chip" flag (done)
# TBD: possibly add coll_distance to sys
# TBD: review coll_error implementation in pt_xfm() [first order only]
# TBD: review and implement the grating xfms as sys params

# TBD: clean up input parameters -- note the way that COL_ERR is *2 in setup now

# 29jan01 -- fixed bug in setup (was a3[3,1] =  cost*sinr)

include	<math.h>
include	<imhdr.h>
include	"deimos.h"
include	"instrument.h"

# define	D_1	20023.15D0
# define	R_CURV	2071.88D0
define	D_1	20018.4D0	# see deimos.h PPLDIST
define	R_CURV	2124.71D0	# ditto; r=83.65in in this case = R_IMSURF


procedure	t_trace()

## char	input[SZ_FNAME]			# input file name (x,y,w)
double	t1, p1
double	t2, p2
double	t3, p3, o3
double	roll3
## pointer	fda

double	e1[3,3], a2[3,3], a3[3,3], a4[3,3]	# transforms
real	ccd[NCCD,3]				# CCD geometry
double	sys[NPARAM]				# system parameters

int	i	# TMP?

int	clgeti()
real	clgetr()
## pointer	open()

begin
	t1 = clgetr ("coll_angle")
	p1 = clgetr ("coll_phi")
	t2 = clgetr ("t2")
	p2 = clgetr ("p2")
	t3 = clgetr ("mu")
	p3 = clgetr ("p3")		# z-misalignment of tilt angle
	o3 = clgetr ("o3")		# independent yaw
	roll3 = clgetr ("roll3")	# roll, assuming p3=0

	CAM_FOC(sys) = clgetr ("cam_foc")
	ORDER(sys) = clgeti ("norder")
	GRLINES(sys) = 1.e-3 * clgeti ("gmm")
	X_OPT(sys) = clgetr ("x_optaxis")
	Y_OPT(sys) = clgetr ("y_optaxis")
	MOS_ROT(sys) = DEGTORAD(clgetr ("mos_rotation"))
# TMP HARDCODES:

COL_DST(sys) = clgetr ("coll_zdst")
COL_ERR(sys) = DEGTORAD(t1)
COL_PHI(sys) = DEGTORAD(p1)

TNT_ANG(sys) = DEGTORAD(71.5 + t2)
TNT_PHI(sys) = DEGTORAD(90. + p2)

# ADDED:
XFO_TLT(sys) = DEGTORAD (0.)
YFO_TLT(sys) = DEGTORAD (0.) - 1.62e-4

MU(sys) = DEGTORAD (t3)
GR_YERR(sys) = DEGTORAD (roll3)
GR_ZERR(sys) = DEGTORAD (o3)
CAM_ANG(sys) = DEGTORAD (clgetr ("cam_angle"))
CAM_PHI(sys) = DEGTORAD (clgetr ("cam_phi"))

do i = 1, 8 {
	CN_XERR(sys,i) = 0.
	CN_YERR(sys,i) = 0.
	CN_RERR(sys,i) = 0.
}



##	call clgstr ("input", input, SZ_FNAME)
##	fda = open (input, READ_ONLY, TEXT_FILE)

	call printf ("# ord=%2f  lines=%5f  tilt=%-6.2f  phi=%-8.4f  roll=%-7.4f yaw=%-7.4f\n")
		call pargd (ORDER(sys))
		call pargd (GRLINES(sys)*1.e3)
		call pargd (t3)
		call pargd (p3)
		call pargd (roll3)
		call pargd (o3)
	call printf ("#  t1=%6.4f p1=%8.4f   t2=%7.4f p2=%8.4f \n")
		call pargd (t1)
		call pargd (p1)
		call pargd (t2)
		call pargd (p2)


	call setup (e1, a2, a3, a4, ccd, sys)


	call trace (e1, a2, a3, a4, ccd, sys)

## Following is the older version
##	call full_trace (fda, e1, a2, a3, a4, ccd, sys)

##	call close (fda)

end

procedure	setup (e1, a2, a3, a4, ccd, sys)

# double	t1, p1		# theta, phi of normal 1
# double	t2, p2		# theta, phi of normal 2
# double	p3, o3	# theta, phi of normal 3
# double	roll3

double	a2[3,3], a3[3,3], a4[3,3]	# transforms
double	e1[3,3]
real	ccd[NCCD,3]	
double	sys[NPARAM]		# system parameters


double	phin, thetan
double	cost, sint, cosp, sinp
double	xsi, cosx, sinx
double	rhon
double	cosr, sinr

begin

# COLL ERROR (first order):
# In this case, we must remove the phi we put in, hence the more complex form
	thetan = 2.*COL_ERR(sys)
	phin = COL_PHI(sys) + HALFPI
	cosp = cos (phin)
	sinp = sin (phin)
	cost = cos (thetan)
	sint = sin (thetan)

	e1[1,1] =  1 - (1-cost)*sinp*sinp	# cosp*cosp + cost*sinp*sinp
	e1[1,2] =  (1-cost) * cosp*sinp
	e1[1,3] = -sinp*sint
	e1[2,1] =  (1-cost) * cosp*sinp
	e1[2,2] =  1 - (1-cost)*cosp*cosp	# sinp*sinp + cost*cosp*cosp
	e1[2,3] =  cosp*sint
	e1[3,1] =  sint*sinp
	e1[3,2] = -sint*cosp
	e1[3,3] =  cost


# TENT MIRROR:
# this mirror is OK to leave in del-theta,phi
	thetan = TNT_ANG(sys)
	phin = TNT_PHI(sys) + HALFPI
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

# MIRROR/GRATING:
# Better in thetax, thetay description (see above)
# for GRATING: need to add the third rotation
# Note the hack with phin + PI.  This is needed to keep the transformed
# x, y axes from flipping, wrt the original x,y.  Not a problem wrt reflections
# but if we want to work _within_ a system it is a problem. Cf. camera also
# note that this _forces_ us to work in a particular hemisphere, and thus we
# must make use of the negative theta as needed.

	thetan = MU(sys)

#	phin = DEGTORAD(p3) + HALFPI
#	if (p3 > 0.) {
#		phin = phin + PI
#		thetan = -thetan
#	}
#
#	cosp = cos (phin)
#	sinp = sin (phin)
#	cost = cos (thetan)
#	sint = sin (thetan)
#	cosx = cos (xsi)
#	sinx = sin (xsi)
#
#	a3[1,1] =  cosx*cosp - cost*sinp*sinx	#  cosp
#	a3[1,2] =  cosx*sinp + cost*cosp*sinx	#  sinp
#	a3[1,3] =  sinx*sint			#  0.
#	a3[2,1] = -sinx*cosp - cost*sinp*cosx	# -cost*sinp
#	a3[2,2] = -sinx*sinp + cost*cosp*cosx	#  cost*cosp
#	a3[2,3] =  cosx*sint			#  sint
#	a3[3,1] =  sint*sinp
#	a3[3,2] = -sint*cosp
#	a3[3,3] =  cost
#
#
#	thetan = -DEGTORAD(t3)
#	thetan = -MU(sys)
#	phin = 0.

#	cosp = 1.
#	sinp = 0.
#	cost = cos (thetan)
#	sint = sin (thetan)
##
#	a3[1,1] =  cosx
#	a3[1,2] =  cost*sinx
#	a3[1,3] =  sinx*sint
#	a3[2,1] = -sinx
#	a3[2,2] =  cost*cosx
#	a3[2,3] =  cosx*sint
#	a3[3,1] =  0.
#	a3[3,2] = -sint
#	a3[3,3] =  cost
#

	thetan = -MU(sys)
	xsi = GR_ZERR(sys)
	rhon = GR_YERR(sys)
	cost = cos (thetan)
	sint = sin (thetan)
	cosx = cos (xsi)
	sinx = sin (xsi)
	cosr = cos (rhon)
	sinr = sin (rhon)

# ... below assumes phin=0. (ie adopts the roll/yaw approach)
	a3[1,1] =  cosx*cosr
	a3[1,2] =  sint*sinr + cost*sinx*cosr
	a3[1,3] = -cost*sinr + sint*sinx*cosr
	a3[2,1] = -sinx
	a3[2,2] =  cost*cosx
	a3[2,3] =  sint*cosx
	a3[3,1] =  cosx*sinr
	a3[3,2] = -sint*cosr + cost*sinx*sinr
	a3[3,3] =  cost*cosr + sint*sinx*sinr

# CAMERA
# again, tranformation from theta-x,y is better, although this is not ridiculous
	thetan = CAM_ANG(sys)
	phin = CAM_PHI(sys) + HALFPI
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

# CCD: define the geometry of the mosaic
	call ccd_geom (ccd, sys)

end


# CCDGEOM: specify the geometry of the CCDs in the mosaic
# Probably want to add FCS devices also (note that FDS CCDs are rotated 90deg)
# Order is x(pix), y(pix), theta(deg)

procedure ccd_geom (a, sys)

real	a[NCCD,3]
double	sys[NPARAM]		# system params

int	i
real	xdimeff, ydimeff	# effective dimensions of CCDs in pixels

begin
	xdimeff = CCDXPIX + 2.*CCDXEDG/PIX_SZ + NOMXGAP/PIX_SZ	# 2135.200
	ydimeff = CCDYPIX + 2.*CCDYEDG/PIX_SZ + NOMYGAP/PIX_SZ	# 4112.000

#	coeff  = pix-off + nom.gap	+ adjustment

#Initial RED Mosaic 02apr24:
#variable plate scale in y, RED Mosaic 02apr26:
# Analysis of first cal sequence, 02apr28
#   1:    0.67     0.34  -0.0033
#   2:    0.52    -0.09  -0.0011
#   3:    0.00     0.00   0.0000
#   4:   -0.07    -0.04  -0.0006
#   5:    0.52    -0.21   0.0030
#   6:    0.39    -0.21   0.0000
#   7:    0.15    -0.24  -0.0043
#   8:   -0.19    -0.52   0.0012

	a[1,1] = -1.5 * xdimeff		+ -20.05
	a[1,2] = -0.5 * ydimeff		+  14.12
	a[1,3] = 0.			+ DEGTORAD(-0.082)

	a[2,1] = -0.5 * xdimeff		+ -12.64
	a[2,2] = -0.5 * ydimeff		+  7.25
	a[2,3] = 0.			+ DEGTORAD(0.030)

	a[3,1] =  0.5 * xdimeff
	a[3,2] = -0.5 * ydimeff
	a[3,3] = 0.

	a[4,1] =  1.5 * xdimeff		+  -1.34
	a[4,2] = -0.5 * ydimeff		+ -19.92
	a[4,3] = 0.			+ DEGTORAD(-0.1206)

	a[5,1] = -1.5 * xdimeff		+ -19.02
	a[5,2] =  0.5 * ydimeff		+  16.46
	a[5,3] = 0.			+ DEGTORAD(0.136)

	a[6,1] = -0.5 * xdimeff		+  -9.65
	a[6,2] =  0.5 * ydimeff		+   8.95
	a[6,3] = 0.			+ DEGTORAD(-0.060)

	a[7,1] =  0.5 * xdimeff		+   1.88
	a[7,2] =  0.5 * ydimeff		+   1.02
	a[7,3] = 0.			+ DEGTORAD(-0.019)

	a[8,1] =  1.5 * xdimeff		+   4.81
	a[8,2] =  0.5 * ydimeff		+ -24.01
	a[8,3] = 0.			+ DEGTORAD(-0.082)

# CN_XOFF(sys,1) = a[1,1]
# CN_YOFF(sys,1) = a[1,2]
# CN_ROT(sys,1)  = a[1,3]

	do i = 1, NCCD {
		a[i,1] = a[i,1] + CN_XERR(sys,i)
		a[i,2] = a[i,2] + CN_YERR(sys,i)
		a[i,3] = a[i,3] + CN_RERR(sys,i)
	}
end


# FULL_TRACE: trace chief rays through system.  The collimator transform must
# be worked out for each ray, so only the collimator error appears as input.
# There will be a simpler version put forth that uses a polynomial mapping
# of input alphas calculated once the errors are established.

procedure	full_trace (fda, e1, a2, a3, a4, ccd, sys)

pointer	fda				# file descriptor for input file
double	e1[3,3]				# collimator error transform
double	a2[3,3], a3[3,3], a4[3,3]	# transforms
real	ccd[NCCD,3]	
double	sys[NPARAM]			# system parameters

real	xpix, ypix

int	n				# CCD number

double	x, y, w

int	stat

int	mosim_coord()
int	fscan(), nscan()

begin

# Ready to loop through:
	while (fscan (fda) != EOF) {

		call gargd (x)
		call gargd (y)
		call gargd (w)

		if (nscan() < 2)
			next

		call pt_xfm (x, y, w, e1, a2, a3, a4, ccd, sys, xpix, ypix, n)
		call ics_to_ccd (xpix, ypix, ccd, n, xpix, ypix)

# Convert to full mosaic image
		stat = mosim_coord (xpix, ypix, n)

		call printf ("%8.2f %8.2f  #  %7.1f \n")
			call pargr (xpix)
			call pargr (ypix)
			call pargd (w*10000.)
	}
end

#
# PT_XFM: point transform; from point on slitmask to point in image.
# Wave < 0 solves for a ghost image.
#

procedure	pt_xfm (xmm, ymm, wave, e1, a2, a3, a4, ccd, sys, xpix, ypix, n)

double	xmm, ymm, wave
double	e1[3,3], a2[3,3], a3[3,3], a4[3,3]	# transforms
real	ccd[NCCD,3]	
double	sys[NPARAM]				# system parameters
real	xpix, ypix
int	n

double	a1[3,3]

double	x, y, pa
double	phi, theta
double	phin, thetan
double	cost, sint, cosp, sinp

double	r[3]

double	alpha, beta, gamma
double	rp, hm
double	tx, ty	# TMP!

real	xp, yp

int	ident_ccd()
begin
	pa = 0.
	x = xmm
	y = ymm
	call mask_to_proj (x, y, pa, x, y, pa)

# comment out above call; add:
# call eprintf ("SPECIAL FOR ZEMAX!!\n")
# refwave=0.5000

# convert to r[3]:
	rp = sqrt (x*x + y*y)
	hm = R_CURV - sqrt (R_CURV*R_CURV - rp*rp)
	theta = atan (rp / (D_1-hm))
	phi = atan2 (y, x)

	cosp = cos (phi)
	sinp = sin (phi)
	cost = cos (theta)
	sint = sin (theta)

	r[1] = cosp * sint
	r[2] = sinp * sint
	r[3] = cost

# call eprintf ("\n[R]: %6f %6f %6f  ")
# call pargd (r[1])
# call pargd (r[2])
# call pargd (r[3])

# COLLIMATOR:
	call coll_angle (x, y, sys, theta, phi)

	thetan = theta
	phin = phi + HALFPI
	cosp = cos (phin)
	sinp = sin (phin)
	cost = cos (thetan)
	sint = sin (thetan)

	a1[1,1] = cosp
	a1[1,2] = sinp
	a1[1,3] = 0.
	a1[2,1] = -cost*sinp
	a1[2,2] = cost*cosp
	a1[2,3] = sint
	a1[3,1] = sint*sinp
	a1[3,2] = -sint*cosp
	a1[3,3] = cost

	call refl (r, a1)
# REVIEW implementation below
	call gen_xfm (r, e1, YES)
	call refl (r, a2)
## NB COMMENT LINE BELOW for ZEMAX cmp
	call gen_xfm (r, a3, YES)

	tx = atan2 (-r[1], -r[3])
	ty = atan2 (-r[2], -r[3])

## HERE IS THE GRATING EQUATION:
## m*wave = -n * space * cos(gamma) * ( sin(beta) + sin(alpha) )

	alpha = -ty
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))

## add the following for ZEMAX compare and comment out remainder of loop
# call printf ("%8.5f %8.5f %8.3f \n")
# call pargd (r[1]/r[3])
# call pargd (r[2]/r[3])
# call pargd (hm)

	beta = asin ((ORDER(sys) * GRLINES(sys) * abs (wave) / cos (gamma)) - sin (alpha))

# If ghost, reflect about grating normal
	if (wave < 0.D0) {
		gamma = -gamma
		beta = -beta
	}

# It is unclear to me right now why gamma does not get changed in sign ...

	call grat_to_ics (beta, gamma, sys, a4, a3, xpix, ypix)

	n = ident_ccd (xpix, ypix)
#	call ics_to_ccd (xpix, ypix, ccd, n, xpix, ypix)
end


# MOSIM_COORD: convert chip coord --> standard mosaic image coordinates
# TMP! TMP!  hardcoded values ...

int	procedure mosim_coord (xpix, ypix, n)

real	xpix, ypix
int	n

int	mx, my
real	flagx, flagy

begin
	if (n > 4) {
		my = 2
		mx = n - 4
	} else {
		my = 1
		mx = n
	}

# flag if off chip:
	if (xpix < 1. || xpix > 2048.)
		flagx = -1.
	else
		flagx = 1.

	if (ypix < 1. || ypix > 4096.)
		flagy = -1.
	else
		flagy = 1.

	xpix = xpix + ((mx-1)*2048 + PRECOL)
	ypix = ypix + (my-1)*4096

	xpix = xpix * flagx
	ypix = ypix * flagy

	if (flagx > 0. && flagy > 0.)
		return (ON_CHIP)
	else
		return (OFF_CHIP)
end


################### FROM NEWLRIS: ################################

define		X_HWID	1.5		# x half-width of triangle function
define		Y_HWID	1.5		# y half-width of triangle function
define		X_CRAD	13		# x centering radius
define		Y_CRAD	13		# y centering radius


# TRACE: This hacked from NEWLRIS version to work with current system.
# Box-centering subroutines borrowed intact.

procedure trace (e1, a2, a3, a4, ccd, sys)

double	e1[3,3], a2[3,3], a3[3,3], a4[3,3]	# transforms
real	ccd[NCCD,3]				# CCD geometry
double	sys[NPARAM]				# system parameters

# new for mosaic images:
pointer	im0			# image pointer to primary HDU
pointer	mmap			# pointer to mosaic vectors
int	nccd

char	slits[SZ_FNAME]			# slit info file
char	output[SZ_FNAME]		# output extraction info
char	image[SZ_FNAME]			# optional image to measure
char	wavelist[SZ_FNAME]		# optional wavelength list
real	refwave				# wavelength
int	index				# identifying index number
real	xsz, ysz			# size in mm of aperture

#
bool	measure				# measure image?
bool	multiwav			# is there a wavelength list?
bool	outfile				# is there an output file?
pointer	im				# image descriptor
pointer	fda, fdb			# file descriptor for in/output files
pointer	fdw				# file descriptor for wavelist

int	nslit				# No. slits
pointer	bufxmm				# pointer to slit x
pointer	bufymm				# pointer to slit y 
pointer	bufwav				# pointer to wavelengths
pointer	bufx, bufy			# pointers to x,y (pixel) slit values
pointer	bufc				# pointer to chip number

char	tchar
char	idstr[SZ_LINE]
int	i, j
int	n
int	nwave
int	nline

# int	niter
# real	xtarg

real	lambda, alpha, beta
real	relscale			# rate
real	anamf				# anamorphic factor in x
# real	xmin, xmax
int	ncol
# real	x1, x2, y1, y2
real	xics, yics		# TMP for FCS work

int	stat
int	mosim_coord()

bool	strne()
int	access(), clgeti()
int	line_count()
int	fscan(), nscan()
real	clgetr()
pointer	open()

begin
# Read in parameters:
	xsz = clgetr ("xpx")
	ysz = clgetr ("ypx")
#	xsz = 4.2
#	ysz = 4.2

	index = clgeti ("index")

	call clgstr ("slits", slits, SZ_FNAME)
	fda = open (slits, READ_ONLY, TEXT_FILE)

	call clgstr ("wavelist", wavelist, SZ_FNAME)
	multiwav = strne (wavelist, "")

#   first count the entries ...
	if (multiwav) {
		fdw = open (wavelist, READ_ONLY, TEXT_FILE)
		nwave = line_count (fdw)
		call malloc (bufwav, nwave, TY_DOUBLE)

		nwave = 0
		while (fscan (fdw) != EOF) {
			call gargwrd (tchar, 1)
			if (tchar == '#') {
				next
			}
			call reset_scan()
			call gargr (refwave)
			Memd[bufwav+nwave] = refwave / 1.e4	# to microns
			if (nscan() < 1)
				next
			nwave = nwave + 1
		}
		call close (fdw)
call eprintf ("Number of input wavelengths: %d\n")
call pargi (nwave)

	} else {
		refwave = clgetr ("refwave")
		nwave = 1
		call malloc (bufwav, nwave, TY_DOUBLE)
		Memd[bufwav] = refwave/1.e4
	}

# Open output file, if requested (FIX -- lots of cleanup ...)
	call clgstr ("output", output, SZ_FNAME)
	outfile = strne (output, "")
	if (outfile && access (output, 0, 0,) == YES) {
#   first count the entries ...
		fdb = open (output, READ_ONLY, TEXT_FILE)
		call close (fdb)
#  ... then reopen
		call eprintf ("appending to %s \n")
			call pargstr (output, SZ_FNAME)
		fdb = open (output, APPEND, TEXT_FILE)
	} else if (outfile) {
		fdb = open (output, NEW_FILE, TEXT_FILE)
	}
	nline = 0

	
# Read in slit positions here
	nslit = line_count (fda)
	call malloc (bufxmm, nslit, TY_DOUBLE)
	call malloc (bufymm, nslit, TY_DOUBLE)
	nslit = 0
	while (fscan (fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#') {
			next
		}
		call reset_scan()
		call gargwrd (idstr, SZ_LINE)
		call gargd (Memd[bufxmm+nslit])
		call gargd (Memd[bufymm+nslit])
		if (nscan() < 3)
			next
		nslit = nslit + 1
	}
	call close (fda)
call eprintf ("Number of input slits: %d\n")
call pargi (nslit)

# Allocate these vectors
	call malloc (bufx, nslit, TY_REAL)	# x (pixels) of slit image
	call malloc (bufy, nslit, TY_REAL)	# y (pixels) of slit image
	call malloc (bufc, nslit, TY_INT)	# CCD number of slit image

call eprintf ("FIX! no alpha/beta! \n")
alpha = DEGTORAD(44.1) - CAM_ANG(sys) + MU(sys)
beta =  CAM_ANG(sys) - MU(sys)

# Solve for scales/anamorphic factor
	if (ORDER(sys) != 0) {
		lambda = 1.e4 / GRLINES(sys) * ORDER(sys) *
			(sin (alpha) - sin (beta))
		relscale = 1.e4 / GRLINES(sys) * ORDER(sys) *
				(cos (alpha) + cos (beta))
		relscale = RADTODEG (1. / relscale)
	} else {
		lambda = INDEF
		relscale = INDEF
	}
		
call eprintf ("\nOn-axis values: beta =%8.4f ==> lambda =%8.2f\n")
call pargr (RADTODEG(beta))
call pargr (lambda)
call eprintf ("Correct at %8.5f degree / A\n")
call pargr (relscale)

	anamf = cos (alpha) / cos (beta)
	call eprintf ("box size:  %4.1f x %4.1f (px)\n")
		call pargr (xsz)
		call pargr (ysz*anamf)

# Measure the centers on an image, if image name given
	call clgstr ("image", image, SZ_FNAME)
	measure = strne (image, "")
	if (!measure) {
		call printf ("No Image for Measured Centers\n")
#####
# TEMPORARY until below is implemented
do i = 0, nslit-1 {
	call pt_xfm (Memd[bufxmm+i], Memd[bufymm+i],
		Memd[bufwav+j], e1, a2, a3, a4, ccd, sys,
		Memr[bufx+i], Memr[bufy+i], n)
	call ics_to_ccd (Memr[bufx+i], Memr[bufy+i], ccd, n, Memr[bufx+i], Memr[bufy+i])
	Memi[bufc+i] = n

## TMP for FCS work:
#		call ccd_to_ics (Memr[bufx+i], Memr[bufy+i], ccd, n, xics, yics)
#		call printf ("IN ICS: %7.3f %7.3f\n")
#			call pargr (xics)
#			call pargr (yics)
## END TMP for FCS work:

# Convert to full mosaic image
	stat = mosim_coord (Memr[bufx+i], Memr[bufy+i], n)
}	
#####
		do i = 0, nslit-1 {
		    call printf ("%03d  %7.3f %7.3f  %7.2f    %7.2f %7.2f \n")
		    	call pargi (index)
		    	call pargd (Memd[bufxmm+i])
		    	call pargd (Memd[bufymm+i])
		    	call pargr (refwave)
		    	call pargr (Memr[bufx+i])
		    	call pargr (Memr[bufy+i])
		}
##		niter = clgeti ("niter")
##		xtarg = clgetr ("xtarget")
##		call xiter_mu (Memd[bufxmm], Memd[bufymm], refwave, xtarg,
##							niter, sys)
	} else {

# Open image:
#		im = immap (image, READ_ONLY, 0)
#		ncol = IM_LEN (im,1)
# Open mosaic image:
		call mos_init (image, "deimos$test/prop.txt", im0, mmap, nccd)

# Loop and print as needed
		do j = 0, nwave-1 {

call eprintf ("processing %d:  %7.2f \n")
call pargi (j)
call pargd (Memd[bufwav+j]*1.D4)

## TMP - FIX THIS
##		    call pt_geom (xmin, 167.8, Memr[bufzxo], Memr[bufzyo],
##			x1, y1, 1, fdz, cd, sys, refwave, xfm)
##		    call pt_geom (xmax, 167.8, Memr[bufzxo], Memr[bufzyo],
##			x2, y2, 1, fdz, cd, sys, refwave, xfm)
##		    if (max (x1, x2) < 1. || min (x1, x2) >= ncol) {
##				call eprintf ("Line off detector\n")
##				next
##		    }

		    do i = 0, nslit-1 {
# added 01/09/05 -- if ghost, check that primary is on chip:
			if (Memd[bufwav+j] < 0.D0) {
			    call pt_xfm (Memd[bufxmm+i], Memd[bufymm+i],
				-1.D0*Memd[bufwav+j], e1, a2, a3, a4, ccd, sys,
				Memr[bufx+i], Memr[bufy+i], n)
			    call ics_to_ccd (Memr[bufx+i], Memr[bufy+i], ccd, n,
						Memr[bufx+i], Memr[bufy+i])
			    stat = mosim_coord (Memr[bufx+i], Memr[bufy+i], n)
			    if (stat == OFF_CHIP) {
				Memr[bufx+i] = INDEF
				next
			    }
			}

			call pt_xfm (Memd[bufxmm+i], Memd[bufymm+i],
				Memd[bufwav+j], e1, a2, a3, a4, ccd, sys,
				Memr[bufx+i], Memr[bufy+i], n)
			call ics_to_ccd (Memr[bufx+i], Memr[bufy+i], ccd, n,
						Memr[bufx+i], Memr[bufy+i])
			Memi[bufc+i] = n

# Convert to full mosaic image
			stat = mosim_coord (Memr[bufx+i], Memr[bufy+i], n)
# added 01/09/04:
			if (stat == OFF_CHIP)
				Memr[bufx+i] = INDEF
		    }	

		    call mbx_recenter (mmap, Memr[bufx], Memr[bufy], Memi[bufc], nslit, xsz, ysz*anamf)
# write measured positions to output file
#	call fprintf (fdb, "# xmm ymm   xpix ypix    ID\n")
		    do i = 0, nslit-1 {
			if (Memr[bufx+i] == INDEF)
				next
			call fprintf (fdb, "%03d %8.3f %8.3f %8.2f    %7.2f %7.2f \n")
				call pargi (index)
				call pargd (Memd[bufxmm+i])
				call pargd (Memd[bufymm+i])
				call pargd (Memd[bufwav+j]*1.D4)
				call pargr (Memr[bufx+i])
				call pargr (Memr[bufy+i])
			 nline = nline + 1
		    }
#		    nline = nline + nslit
		}

# Write the setup info:
	call fprintf (fdb, "# %03d %s %7.3f  %6.3f %6.3f  %5.0f %3.0f %4d\n")
		call pargi (index)
		call pargstr (image)
		call pargd (RADTODEG(MU(sys)))
		call pargd (RADTODEG(GR_YERR(sys)))
		call pargd (RADTODEG(GR_ZERR(sys)))
		call pargd (GRLINES(sys)*1.e3)
		call pargd (ORDER(sys))
		call pargi (nline)


# Close up
#		call imunmap (im)
# TMP! FIX! Need a mos_unmap
	}


# PROOF section
if (multiwav && !measure && outfile) {
	call proof (Memd[bufxmm], Memd[bufymm], nslit, Memd[bufwav], nwave, fdb,
			e1, a2, a3, a4, ccd, sys)
}


	if (outfile)
		call close (fdb)

# Release memory
	call mfree (bufc, TY_INT)
	call mfree (bufy, TY_REAL)
	call mfree (bufx, TY_REAL)
	call mfree (bufymm, TY_DOUBLE)
	call mfree (bufxmm, TY_DOUBLE)
	call mfree (bufwav, TY_DOUBLE)
end


#
# MBX_RECENTER:
# Changes from previous:
# 0. mmap instead of im pointer, chip vector
# 1. added mosaic/chip_section()

procedure	mbx_recenter (mmap, xs, ys, cs, nslit, xsz, ysz)

pointer	mmap				# mosaic image structure
real	xs[nslit], ys[nslit]		# x,y vectors of slit/boxes
int	cs[nslit]			# CCD number
int	nslit
real	xsz, ysz			# Full size of boxes (pix)

int	nxbox, nybox			# Size of subraster

int	i, j, k

int	nx, ny
real	xb, yb
int 	x1, y1, x2, y2
pointer	bufx, bufy, bufzx, bufzy
pointer	buf

pointer	chip_sect()

int	clgeti()
begin
	nxbox = clgeti ("nxbox")
	nybox = clgeti ("nybox")

	nx = nxbox		# fixed length now that chip_sect zero-fills
	ny = nybox

# Allocate arrays for marginal plots
	call malloc (bufx, nxbox, TY_REAL)
	call malloc (bufy, nybox, TY_REAL)
	call malloc (bufzx, nxbox, TY_REAL)
	call malloc (bufzy, nybox, TY_REAL)

# Find the Box; we do a kludge -- run this twice, recentered on box the 2nd time
# probably not necessary
	do k = 1, nslit {
		xb = xs[k]
		yb = ys[k]
# mod. 01/09/04:
#		if (xb == INDEF || xb < 1.+xsz || xb > ncols-xsz ||
#					yb < 1.+ysz || yb > nlines-ysz) {
		if (xb == INDEF) {
			xs[k] = INDEF
			ys[k] = INDEF
			next
		}

		do j = 1, 2 {
		    x1 = xb - nxbox/2
		    x2 = x1 + nxbox - 1
		    y1 = yb - nybox/2
		    y2 = y1 + nybox - 1

# Get the image section (limit checking and zero-fill included)
		    buf = chip_sect (x1, x2, y1, y2, mmap, cs[k], YES)

# Fill position vectors
		    do i = 0, nx-1 {
			Memr[bufx+i] = i + x1
		    }
		    do i = 0, ny-1 {
			Memr[bufy+i] = i + y1
		    }

		    call box_center (Memr[bufx], Memr[bufy], Memr[buf], nx, ny,
			Memr[bufzx], Memr[bufzy], xsz, ysz, xb, yb)

		    if (xb == INDEF || yb == INDEF)
			break
		}

		call printf ("Box center:  %6.2f %6.2f  (del:%4.1f,%4.1f)\n")
			call pargr (xb)
			call pargr (yb)
			if (xb != INDEF && yb != INDEF) {
				call pargr (xb-xs[k])
				call pargr (yb-ys[k])

# added 01/09/05 TMP!! HARDCODE!!
			    if ((xb-xs[k])**2 + (yb-ys[k])**2 > 100.)
				xb = INDEF

			} else {
				call pargr (INDEF)
				call pargr (INDEF)
			}

		xs[k] = xb		# replace with new coord.
		ys[k] = yb		# replace with new coord.
	}

# Done with the box-finding: clean up
	call mfree (bufzy, TY_REAL)
	call mfree (bufzx, TY_REAL)
	call mfree (bufy, TY_REAL)
	call mfree (bufx, TY_REAL)
end

procedure	box_center (x, y, z, nx, ny, zx, zy, xsz, ysz, cx, cy)

real	x[nx], y[ny]		# position vectors
real	z[nx,ny]		# intensity array
int	nx, ny			# size of vectors/array
real	zx[nx], zy[ny]		# Work vectors for x,y cuts
real	xsz, ysz		# pixel sizes of boxes
real	cx, cy			# returned centers
real	grad			# new gradient, currently unused

int	i, j, i1, i2, j1, j2

real	saw_xcorr()
begin

# Get the initial profile
	call amovkr (0., zx, nx)
	call amovkr (0., zy, ny)
	j1 = 0.5 * ny - 1
	j2 = j1 + 4
	do j = j1, j2 {
		call aaddr (zx, z[1,j], zx, nx)
	}
	i1 = 0.5 * nx - 1
	i2 = i1 + 4
	do j = 1, ny {
	    do i = i1, i2 {
		zy[j] = zy[j] + z[i,j]
	    }
	}
	call amulkr (zx, 1./real(j2-j1+1), zx, nx)
	call amulkr (zy, 1./real(i2-i1+1), zy, ny)

	cx = saw_xcorr (x, zx, nx, xsz, X_HWID, X_CRAD, grad)
	if (cx ==INDEF) {
		call eprintf ("Can't find box at %4.0f,%4.0f\n")
			call pargr (x[nx/2])
			call pargr (y[ny/2])
	}
	cy = saw_xcorr (y, zy, ny, ysz, Y_HWID, Y_CRAD, grad)
	if (cy ==INDEF) {
		call eprintf ("Can't find box at %4.0f,%4.0f -- check positions on image and prescan value!\n")
			call pargr (x[nx/2])
			call pargr (y[ny/2])
	}

	if (cx == INDEF || cy == INDEF)
		return

# Get the final profile:
	call amovkr (0., zx, nx)
	call amovkr (0., zy, ny)
	j1 = cy - y[1] + 1 - 0.5 * xsz + X_HWID
	j2 = cy - y[1] + 1 + 0.5 * xsz - X_HWID + 0.5
	do j = j1, j2 {
		call aaddr (zx, z[1,j], zx, nx)
	}
	i1 = cx - x[1] + 1 - 0.5 * xsz + X_HWID
	i2 = cx - x[1] + 1 + 0.5 * xsz - X_HWID + 0.5
	do j = 1, ny {
	    do i = i1, i2 {
		zy[j] = zy[j] + z[i,j]
	    }
	}
	call amulkr (zx, 1./real(j2-j1+1), zx, nx)
	call amulkr (zy, 1./real(i2-i1+1), zy, ny)

	cx = saw_xcorr (x, zx, nx, xsz, X_HWID, X_CRAD, grad)
	if (cx ==INDEF)
		call eprintf ("Got lost on recenter; check input positions")
	cy = saw_xcorr (y, zy, ny, ysz, Y_HWID, Y_CRAD, grad)
	if (cy ==INDEF)
		call eprintf ("Got lost on recenter; check input positions")

	if (cx == INDEF || cy == INDEF)
		return

# Print out centers for debugging
# call eprintf ("cx,cy final = %6.2f %6.2f\n")
# call pargr (cx)
# call pargr (cy)

end

#
# SAW_XCORR: xcorr's a modified sawtooth with a vector to find centers
# This subroutine ALSO appears in MBOXFIND
#
# NOTE -- there can/will be an error if "sz-2*hwid" is larger than feature.
#
# Note -- as this stands, it is still not quite correct. If there are multiple
# crossings, the check on the slope SHOULD BE over the width of the feature;
# instead it is local only. This will probably produce the desired result in
# realistic cases.

real	procedure saw_xcorr (xgrid, z, nx, sz, hwid, noff, grad)

real	xgrid[nx], z[nx]		# position and intensity vectors
int	nx				# size of above
real	sz				# spacing of the peaks
real	hwid				# half-width of the peaks (pixels)
int	noff				# pixels to correlate across
real	grad				# returned gradient

int	nf, nc, xcen
int	ipos, ineg
int	i
real	xpeak, x, xdiff
real	maxderiv, cpos, cneg, fractpos
pointer	form, xcorr, deriv

real	adotr()

begin
# Create "form" vector; this is offset pos and neg. triangle functions
#          +
#         + +
#  +++++++   +++++++++   +++++++++.
#                     + +
#                      +
#
	nf = nx + 2 * noff
	nc = 1  + 2 * noff
	call calloc (form, nf, TY_REAL)
	call calloc (xcorr, nc, TY_REAL)
	call calloc (deriv, nc, TY_REAL)

	xcen = nf / 2 			# Center (to low ) of vectors
	xpeak = 0.5 * sz
	do i = 0, nf-1 {
#		x = i + 1 - (xcen + noff) + 1	# relative offset
		x = nf - (i+1) - xcen		# relative offset
		xdiff = abs (abs (x) - xpeak)
		if (x > 0.)
			Memr[form+i] =  max ((1. - xdiff / hwid), 0.)
		else
			Memr[form+i] = -max ((1. - xdiff / hwid), 0.)
	}

# ready to x-corr; this is actually backwards, but since "form" is symmetric we
# may cheat and assign the negative value:

	do i = 0, 2*noff
		Memr[xcorr+i] = adotr (z, Memr[form+2*noff-i], nx)
	do i = 0, 2*noff-1
		Memr[deriv+i] = 0.5 * (Memr[xcorr+i+1] - Memr[xcorr+i-1])
	Memr[deriv+2*noff] =  Memr[xcorr+2*noff] - Memr[xcorr+2*noff-1]
	Memr[deriv]        =  Memr[xcorr+1]      - Memr[xcorr]

	maxderiv = 0.
	cneg = INDEF
	do i = 0, 2*noff-1 {
		if (Memr[xcorr+i] <= 0. && Memr[xcorr+i+1] > 0.) {
			if (Memr[deriv+i] > maxderiv) {
				ineg = noff - i
				ipos = ineg - 1
				cneg = Memr[xcorr+i]
				cpos = Memr[xcorr+i+1]
				maxderiv = Memr[deriv+i]
			}
		}
	}
	if (cneg == INDEF) {
		x = INDEF
		return (x)
	}

	fractpos = cpos / (cpos - cneg)
	xcen = (nx+1) /2
	x = fractpos * xgrid[xcen-ineg] + (1. - fractpos) * xgrid[xcen-ipos]

	call mfree (deriv, TY_REAL)
	call mfree (xcorr, TY_REAL)
	call mfree (form, TY_REAL)

# TMP! (trial)
	grad = maxderiv

	return (x)	
end



procedure	proof (xmm, ymm, nslit, wave, nwave, fdb, e1, a2, a3, a4, ccd, sys)

double	xmm[nslit], ymm[nslit]
int	nslit
double	wave[nwave]
int	nwave
pointer	fdb				# pointer to output file
double	e1[3,3], a2[3,3], a3[3,3], a4[3,3]	# transforms
real	ccd[NCCD,3]	
double	sys[NPARAM]				# system parameters


char	tchar
pointer	fda
pointer	bufx, bufy			# pointers
pointer bufx0,bufy0
int	ndx
int	i, j
int	n, stat
real	xpix, ypix

real	theta, cost, sint
real	delx, dely, resx, resy
real	xp, yp, xoff, yoff
real	x0cen, y0cen, xcen, ycen, xc, yc

int	mosim_coord()
pointer	open()
int	fscan(), nscan()

begin

# clean out the CCD geom
do i = 1, 8 {
	CN_XERR(sys,i) = 0.
	CN_YERR(sys,i) = 0.
	CN_RERR(sys,i) = 0.
}

# Allocate these vectors
	call malloc (bufx, nslit*nwave, TY_REAL)	# x (pixels) of slit
	call malloc (bufy, nslit*nwave, TY_REAL)	# y (pixels) of slit

	ndx = 0

	do j = 1, nwave {
	    do i = 1, nslit {
		call pt_xfm (xmm[i], ymm[i], wave[j], e1, a2, a3, a4, ccd, sys, xpix, ypix, n)
		call ics_to_ccd (xpix, ypix, ccd, n, xpix, ypix)

# reconvert to ICS
call eprintf ("x,y,n: %5f %5f %d\n")
call pargr (xpix)
call pargr (ypix)
call pargi (n)
		call ccd_to_ics (abs(xpix), abs(ypix), ccd, n,
						Memr[bufx+ndx], Memr[bufy+ndx])

#		stat = mosim_coord (Memr[bufx+ndx], Memr[bufy+ndx], n)
		ndx = ndx + 1
	    }	
	}

	call malloc (bufx0, nslit*nwave, TY_REAL)	# x (pixels) of slit
	call malloc (bufy0, nslit*nwave, TY_REAL)	# y (pixels) of slit

	fda = open ("baseline", READ_ONLY, TEXT_FILE)

call eprintf ("%d * %d \n")
call pargi (nslit)
call pargi (nwave)
	ndx = 0
	while (fscan (fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#') {
			next
		}
		call reset_scan()
		call gargr (Memr[bufx0+ndx])
		call gargr (Memr[bufy0+ndx])
		if (nscan() < 2)
			next
		ndx = ndx + 1
		if (ndx > nslit*nwave)
			call fatal (0, "list problem")
	}

	if (ndx != nslit*nwave)
		call fatal (0, "lists dont match")

# Work out transform (based on specific geometry of first two spots)
	x0cen = 0.5 * (Memr[bufx0] + Memr[bufx0+1])
	y0cen = 0.5 * (Memr[bufy0] + Memr[bufy0+1])
	xcen = 0.5 * (Memr[bufx] + Memr[bufx+1])
	ycen = 0.5 * (Memr[bufy] + Memr[bufy+1])

	xoff = (Memr[bufx] - Memr[bufx0])
	yoff = (Memr[bufy] - Memr[bufy0])
	xoff = 0.5 * (xoff + Memr[bufx+1] - Memr[bufx0+1])
	yoff = 0.5 * (yoff + Memr[bufy+1] - Memr[bufy0+1])

	delx = Memr[bufx+1] - Memr[bufx]
	dely = Memr[bufy+1] - Memr[bufy]
	theta = atan2 (dely, delx)

	delx = Memr[bufx0+1] - Memr[bufx0]
	dely = Memr[bufy0+1] - Memr[bufy0]
	theta = theta - atan2 (dely, delx)

	cost = cos (theta)
	sint = sin (theta)
	
call eprintf ("%5f %5f %5f\n")
call pargr (xoff)
call pargr (yoff)
call pargr (RADTODEG(theta))

	do i = 0, ndx-1 {
		delx = Memr[bufx+i] - Memr[bufx0+i] - xoff
		dely = Memr[bufy+i] - Memr[bufy0+i] - yoff

		xc = Memr[bufx0+i] - x0cen
		yc = Memr[bufy0+i] - y0cen
		xp =  cost * xc - sint * yc + xcen
		yp =  sint * xc + cost * yc + ycen
		resx = Memr[bufx+i] - xp 
		resy = Memr[bufy+i] - yp 
#		call fprintf (fdb, "%8.2f %8.2f  %5.2f %5.2f %6.2f %6.2f \n")
		call fprintf (fdb, "vector %7.1f %7.1f  %6.1f %6.1f \n")
			call pargr (Memr[bufx+i])
			call pargr (Memr[bufy+i])
			call pargr (resx * 250.)
			call pargr (resy * 250.)
#			call pargr (delx*1.)
#			call pargr (dely*1.)
	}

	do i = 0, 1 {
		call fprintf (fdb, "vector %7.1f %7.1f  %6.1f %6.1f \n")
			call pargr (Memr[bufx+i])
			call pargr (Memr[bufy+i]-5000.)
			call pargr ((Memr[bufx0+i] - Memr[bufx+i]) * 25.)
			call pargr ((Memr[bufy0+i] - Memr[bufy+i]) * 25.)
	}

	call mfree (bufx0, TY_REAL)
	call mfree (bufy0, TY_REAL)
	call mfree (bufx, TY_REAL)
	call mfree (bufy, TY_REAL)
end

