# 1-D spectrum; based on TRACE

include <imhdr.h>
include <math.h>
include <math/curfit.h>
include <math/iminterp.h>
include "instrument.h"

define	D_1	20018.4D0	# see deimos.h PPLDIST
define	R_CURV	2124.71D0	# ditto; r=83.65in in this case = R_IMSURF

procedure t_qextr ()

char	image[SZ_FNAME]			# input image
char	outspec[SZ_FNAME]		# output spectrum
double	xslit, yslit			# x,y of slit (mm)
real	objhwid, skip, skywin		# extraction window params
pointer	im, im2

double	e1[3,3], a2[3,3], a3[3,3], a4[3,3]	# transforms
real	ccd[NCCD,3]				# CCD geometry
double	sys[NPARAM]				# system parameters

double	r[3]				# direction vector
double	alpha, beta, gamma		# grating angles
double	sign

real	wave[8192]			# wavelength at each pixel
real	obj[8192], skysum[8192]		# summed obj, sky at each pixel
real	xpix[8192], ypix[8192]		# x,y values at each extraction point
real	win[6]				# window params

double	lam0, lamb, lamr
double	lamgb, lamgr
#double	x, y
double	w, winc
int	nres
real	xccd, yccd			# Input CCD values
real	xics, yics			# returned ICS values
pointer	bufx, bufy			# pointers to trace values
pointer	bufwt				# weights vector for curfit (unused)
pointer	cv1, cv2			# pointers to curfits for low,upp CCDs
pointer	buflin

real	w0, dw			# wavelength scale for output spectrum
int	npix			# number of pix in initial spectrum
int	nx			# number of pix in output spectrum

int	sub_sky				# subtract the sky?
int	i, j
int	line
int	n
#int	nlow
int	n1, n2
int	stat, wt
int	xerr

int	col_extract()
int	mosim_coord()

bool	clgetb()
# bool	strne()
# int	access() 
# int	fscan(), nscan()
int	clgeti()
real	clgetr()
real	cveval()
pointer	immap(), imgs2r(), impl2r()

begin

# Read in parameters:
	call clgstr ("image", image, SZ_FNAME)
	call clgstr ("outspec", outspec, SZ_FNAME)

	xslit = clgetr ("xslit")
	yslit = clgetr ("yslit")

	dw = clgetr ("dw") * 1.e-4	## TMP! FIX -- all in microns

	if (clgetb ("ghost"))
		sign = -1.d0
	else
		sign = 1.D0

COL_DST(sys) = 2197.1
	COL_ERR(sys) = DEGTORAD(clgetr ("trace.coll_angle"))
	COL_PHI(sys) = DEGTORAD(clgetr ("trace.coll_phi"))

TNT_ANG(sys) = DEGTORAD(71.5 + clgetr ("trace.t2"))
TNT_PHI(sys) = DEGTORAD(90. + clgetr ("trace.p2"))

	CAM_ANG(sys) = DEGTORAD(clgetr ("trace.cam_angle"))
	CAM_PHI(sys) = DEGTORAD(clgetr ("trace.cam_phi"))

	CAM_FOC(sys) = clgetr ("trace.cam_foc")

	X_OPT(sys) = clgetr ("trace.x_optaxis")
	Y_OPT(sys) = clgetr ("trace.y_optaxis")
	MOS_ROT(sys) = DEGTORAD(clgetr ("trace.mos_rotation"))

	MU(sys) = DEGTORAD(clgetr ("mu"))
	GR_YERR(sys) = DEGTORAD(clgetr ("roll3"))	# roll, assuming p3=0
	GR_ZERR(sys) = DEGTORAD(clgetr ("o3"))		# independent yaw
	ORDER(sys) = clgeti ("trace.norder")	# grating order
	GRLINES(sys) = 1.e-3 * clgetr ("trace.gmm")	# grating rullings

# Read in window params  ## MUST FIX
	objhwid = 0.5 * clgetr ("objwin")
	skip = clgetr ("skip")
	skywin = clgetr ("skywin")

	win[3] = -objhwid
	win[2] = win[3] - skip
	win[1] = win[2] - skywin
	win[4] = objhwid
	win[5] = win[4] + skip
	win[6] = win[5] + skywin

	if (skywin > 0.)
		sub_sky = YES
	else
		sub_sky = NO

if (sub_sky == YES)
call eprintf ("subtracting sky!\n")


# Read to start: setup and define alpha
	call setup (e1, a2, a3, a4, ccd, sys)

# call debug_mod (e1, a2, a3, a4, ccd, sys)

	call pre_grating (xslit, yslit, e1, a2, sys, r)

# transform into the grating system
	call gen_xfm (r, a3, YES)

	alpha = -atan2 (-r[2], -r[3])
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))

# what is on-axis, are reasonable extremes? (Note 4400 includes gaps, X,Y_OPT)
	xics = 0.
	yics = 0.
	call ics_to_grat (xics, yics, sys, a4, a3, beta, gamma)
	lam0 = 1. / (ORDER(sys) * GRLINES(sys)) * (sin (alpha) + sin (beta))

	xics = -4000.
	yics = -4400.
	call ics_to_grat (xics, yics, sys, a4, a3, beta, gamma)
	lamb = cos (gamma) / (ORDER(sys) * GRLINES(sys)) * (sin (alpha) + sin (beta))

	xics = 0.
	yics = 4400.
	call ics_to_grat (xics, yics, sys, a4, a3, beta, gamma)
	lamr = cos (0.) / (ORDER(sys) * GRLINES(sys)) * (sin (alpha) + sin (beta))

call eprintf ("Wavelength range %5f -- %5f (cen= %5f)\n")
call pargd (lamb)
call pargd (lamr)
call pargd (lam0)

	xics =  0.
	yics = -4400.
	call ics_to_grat (xics, yics, sys, a4, a3, beta, gamma)
	lamgb = cos (0.) / (ORDER(sys) * GRLINES(sys)) * (sin (alpha) + sin (-beta))

	xics = -4400.
	yics = 4400.
	call ics_to_grat (xics, yics, sys, a4, a3, beta, gamma)
	lamgr = cos (gamma) / (ORDER(sys) * GRLINES(sys)) * (sin (alpha) + sin (-beta))

call eprintf ("Ghost Wavelength range %5f -- %5f \n")
call pargd (lamgb)
call pargd (lamgr)

if (sign < 0) {
lamb = lamgb
lamr = lamgr
}

########## BEGIN blaze test  #####
#   call blaze_mod (alpha, sys, lamb, lamr)
########## END blaze test  #####

	call cvinit (cv1, SPLINE3, 11, -10., 4106.)
	call cvinit (cv2, SPLINE3, 11, -10., 4106.)

	nres = 1000
	call malloc (bufx, nres, TY_REAL)
	call malloc (bufy, nres, TY_REAL)
	call malloc (bufwt, nres, TY_REAL)

	winc = (lamr - lamb) / (nres - 1)
	do i = 0, nres-1 {
		w = lamb + i * winc
		w = sign * w
		call pt_xfm (xslit, yslit, w, e1, a2, a3, a4, ccd, sys, xccd, yccd, n)
		Memr[bufx+i] = xccd
		Memr[bufy+i] = yccd

		if (i == 0) {
			n1 = n
		}
		if (n == n1) {
#			nlow = i
			if (yccd < -10. || yccd > 4106.)
				next
			call cvaccum (cv1, yccd, xccd, wt, WTS_UNIFORM)
		} else {
			n2 = n
			call cvaccum (cv2, yccd, xccd, wt, WTS_UNIFORM)
		}
	}


## Fit trace in two parts (FIX Hardcoded order)
	call cvsolve (cv1, stat)
	call cvsolve (cv2, stat)


# Ready to extract spectrum
	im = immap (image, READ_ONLY, 0)
	im2 = immap (outspec, NEW_COPY, im)

	xerr = 0

	do j = 1, 4096 {
		xccd = cveval (cv1, real (j))
		yccd = j
		call ccd_to_grating (xccd, yccd, n1, sys, ccd, a4, r)
		xpix[j] = xccd
		ypix[j] = yccd

# transform backwards into grating system
		call gen_xfm (r, a3, YES)

		beta = -atan2 (r[2], r[3])
		gamma = asin (r[1])
# if (j-(j/20)*20 == 0) {
# call eprintf ("%d= %7f\n")
# call pargi (j)
# call pargd (gamma)
# }

		wave[j] = cos (gamma) / (ORDER(sys) * GRLINES(sys)) * (sin (alpha) + sin (sign*beta))
		stat = mosim_coord (xccd, yccd, n1)
		line = yccd + 0.5
		buflin = imgs2r (im, 1, 8240, line, line)
		xerr = xerr + col_extract (Memr[buflin], 8240, xccd, win,
					sub_sky, obj[j], skysum[j])

#call printf ("%7.2f %8.1f # %6.1f %d # %8.1f %8.1f\n")
#call pargr (wave[j]*1.e4)
#call pargr (obj[j]-skysum[j])
#call pargr (xccd)
#call pargi (j)
#call pargr (obj[j])
#call pargr (skysum[j])
	}

	do j = 4097, 8192 {
		xccd = cveval (cv2, real (j-4096))
		yccd = j - 4096
		call ccd_to_grating (xccd, yccd, n2, sys, ccd, a4, r)
		xpix[j] = xccd
		ypix[j] = j

## transform backwards into grating system
		call gen_xfm (r, a3, YES)

		beta = -atan2 (r[2], r[3])
		gamma = asin (r[1])

		wave[j] = cos (gamma) / (ORDER(sys) * GRLINES(sys)) * (sin (alpha) + sin (sign*beta))
		stat = mosim_coord (xccd, yccd, n2)
		line = yccd + 0.5
		buflin = imgs2r (im, 1, 8240, line, line)
		xerr = xerr + col_extract (Memr[buflin], 8240, xccd, win,
					sub_sky, obj[j], skysum[j])

# call printf ("%7.2f %8.1f # %6.1f %d # %8.1f %8.1f\n")
# call pargr (wave[j]*1.e4)
# call pargr (obj[j]-skysum[j])
# call pargr (xccd)
# call pargi (j)
# call pargr (obj[j])
# call pargr (skysum[j])
	}

	npix = 8192			# NB: note number of hard-codes above!
	if (sub_sky == NO)
		call amovkr (0., skysum, npix)

	w0 = wave[1]
	w0 = w0 - dw * (w0/dw - int (w0/dw))
	nx = (wave[npix] - w0) / dw + 1

	IM_LEN(im2,1) = nx
	IM_LEN(im2,2) = 5
	IM_PIXTYPE(im2) = TY_REAL
# Add relevant keywords to output header
	call imastr (im2, "CTYPE1", "LINEAR")
	call imaddr (im2, "CRPIX1", 1.)
	call imaddr (im2, "CRVAL1", w0*1.e4)
	call imaddr (im2, "CD1_1", dw*1.e4)
	call imaddr (im2, "CDELT1", dw*1.e4)

# Axis 2 wcs should be OK, add axis 3
#	call imastr (im2, "WCSDIM", 3)
#	call imastr (im2, "CTYPE3", "LINEAR")
#	call imaddr (im2, "CRPIX3", 1.)
#	call imaddr (im2, "CRVAL3", 1.)
#	call imaddr (im2, "CDELT3", 1.)
#	call imaddr (im2, "CD3_3", 1.)
#	call imaddr (im2, "LTM3_3", 1.)

# Change system for splot
	call imastr (im2, "WAT0_001", "system=equispec")

call eprintf ("no flux conservation applied\n")

	call lin_spec (wave, obj, npix, Memr[impl2r (im2, 1)], w0, dw, nx, NO)
	call lin_spec (wave, skysum, npix, Memr[impl2r (im2, 2)], w0, dw, nx, NO)
	call asubr (obj, skysum, obj, npix)
	call lin_spec (wave, obj, npix, Memr[impl2r (im2, 3)], w0, dw, nx, NO)
	call lin_spec (wave, xpix, npix, Memr[impl2r (im2, 4)], w0, dw, nx, NO)
	call lin_spec (wave, ypix, npix, Memr[impl2r (im2, 5)], w0, dw, nx, NO)

	xerr = min (xerr, 1)

	call imunmap (im2)
	call imunmap (im)
	call cvfree (cv2)
	call cvfree (cv1)
	call mfree (bufwt, TY_REAL)
	call mfree (bufy, TY_REAL)
	call mfree (bufx, TY_REAL)

end





#
# PRE_GRATING: calculate the direction vector _before_ the grating
#

procedure	pre_grating (xmm, ymm, e1, a2, sys, r)

double	xmm, ymm
double	e1[3,3], a2[3,3]			# pre-grating transforms
double	sys[NPARAM]				# system parameters
double	r[3]

double	a1[3,3]

double	x, y, pa
double	phi, theta
double	phin, thetan
double	cost, sint, cosp, sinp

double	rp, hm


begin
	pa = 0.
	x = xmm
	y = ymm
	call mask_to_proj (x, y, pa, x, y, pa)

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
	call gen_xfm (r, e1, YES)
	call refl (r, a2)
end



#
# CCD_TO_GRATING: traces x, y back to direction vector post-grating.

procedure	ccd_to_grating (xccd, yccd, n, sys, ccd, a, r)

real	xccd, yccd
int	n				# CCD number
double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]				# CCD geometry
double	a[3,3]				# camera tranformation
double	r[3]			# coord vector


double	theta, phi
double	xp, yp, rp
double	xpp, ypp
double	cosp, sinp, cost, sint

real	x, y

begin
	call ccd_to_ics (xccd, yccd, ccd, n, x, y)

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

end

procedure	lin_spec (w, z, nw, zout, w0, dw, nx, flux)

real	w[nw]			# wavelength values
real	z[nw]			# intensity values
int	nw			# number in vector
real	zout[nx]		# vector of linearized output values
real	w0, dw			# starting, specing of output wavelengths
int	nx			# number of pixels in output spectrum
int	flux			# preserve flux?

int	i
pointer	bufwout				# vector of desired wavelengths
pointer	bufwx				# vector of desired wavelength x-vals
pointer	bufnorm				# vector to flux normalization factors
pointer	asi

begin

# Allocate vectors
	call malloc (bufwout, nx, TY_REAL)
	call malloc (bufwx,   nx, TY_REAL)
	call malloc (bufnorm, nx, TY_REAL)

# Initialize the interpolation routines
	call asiinit (asi, II_SPLINE3)		# SINC slow, poor around sat.
		# SPLINE3 shows no sign. improv., but better noise char.

# assign the output wavelength scale
	do i = 0, nx-1
		Memr[bufwout+i] = w0 + i * dw


# Interpolate the position of the output wavelength scale
	call interp_wav (w, nw, Memr[bufwout], Memr[bufwx], nx)

# fit and interpolate input line
	call asifit (asi, z, nw)

	call asivector (asi, Memr[bufwx], zout, nx)

	if (flux == YES) {
# Approx. Flux conservation; flag pixels with _NO_ flux:
		Memr[bufnorm] = 0.
		Memr[bufnorm+nx-1] = 0.
		call amovr (Memr[bufwx+0], Memr[bufnorm+1], nx-2)
		call asubr (Memr[bufwx+2], Memr[bufnorm+1],
							Memr[bufnorm+1], nx-2)
		call amulkr (Memr[bufnorm], 0.5, Memr[bufnorm], nx)

		do i = 0, nx-3 {
			if (Memr[bufnorm+i] != 0.) {
				Memr[bufnorm+i] = Memr[bufnorm+i+2]
				Memr[bufnorm+i+1] = Memr[bufnorm+i+2]
				break
			}
		}
		do i = nx-1, 2, -1 {
			if (Memr[bufnorm+i] != 0.) {
				Memr[bufnorm+i] = Memr[bufnorm+i-2]
				Memr[bufnorm+i-1] = Memr[bufnorm+i-2]
				break
			}
		}

# ... bufnorm is the inv. spacing in A; here convert to units of wavinc
	    call amulkr (Memr[bufnorm], dw, Memr[bufnorm], nx)

	    call amulr (zout, Memr[bufnorm], zout, nx)

	    do i = 0, nx-1 {
		if (Memr[bufnorm+i] == 0.)
			zout[i+1] = 0.
	    }
	}
	

	call asifree (asi)

	call mfree (bufnorm, TY_REAL)
	call mfree (bufwx, TY_REAL)
	call mfree (bufwout, TY_REAL)
end


######################### from REDUX.SHAPEUP #######################

#
# INTERP_WAV
# This procedure estimates the x-position of a vector of desired wavelengths,
# based on a vector of known (or calculated) wavelengths.  The input and output
# vector need NOT be of equal length, and the output grid can overwrite the
# input grid.
# Both grids must be SORTED low-to-high!

procedure interp_wav (wave, nx, wgrid, x, nxout)

real	wave[nx]			# input wavelengths (sorted)
int	nx				# length of input vector
real	wgrid[nxout]			# wavelengths desired (input;sorted)
real	x[nxout]			# these are the x-positions (returned)
int	nxout

int	i, j, jstart
real	w, wmin, wmax

begin
	wmin = wave[1]
	wmax = wave[nx]
	jstart = 1

	do i = 1, nxout {
# check to see if within data limits; set to limits if so
		w = wgrid[i]
		if (w <= wmin) {
			x[i] = 1
			next
		}
		if (w >= wmax) {
			x[i] = nx
			next
		}

# search for correct value in between
		do j = jstart, nx-1 {
			if (w < wave[j+1]) {
			    x[i] = j + (w - wave[j]) / (wave[j+1] - wave[j])
			    jstart = j
			    break
			}
		}
	}
end

######################### from REDUX.SHIPOUT #######################

#
# COL_EXTRACT: get object, scaled-sky counts along one column; does limit checks
# return error on window truncation

int	procedure col_extract (slit, nx, xcen, win, extr_sky, objval, skyval)

real	slit[nx]			# center of window (vector)
int	nx				# length of slit vector
real	xcen				# center of windows
real	win[6]				# window parameters
int	extr_sky			# extract and scale the sky?
real	objval, skyval			# object, sky values (returned)

real	x1, x2, x3, x4, x5, x6		# edges to windows
real	xmaxf, xcenf			# max, cent x in fractional px conv.
int	i, i1, i2
int	xerr				# are windows truncated?

begin

### FRACTIONAL PIXEL CONVENTION:
# We wish to work in pixels that run from 0.0-1.0.  This means pixels specified
# by the usual convention correspond to 0.5px higher value, as adjusted below.
# Also, then, the entire pixel range is [1.0,nx+1.0], with pixel centers
# at 1.5 and nx+0.5 at the extremes.

	xmaxf = nx + 1.0
	xcenf = xcen + 0.5		# 0.5 pix offset for pix-cen convention
	x1 = xcenf + win[1]
	x2 = xcenf + win[2]
	x3 = xcenf + win[3]
	x4 = xcenf + win[4]
	x5 = xcenf + win[5]
	x6 = xcenf + win[6]

	skyval = 0.
	objval= 0.
	xerr = 0
	if (x3 <= xmaxf && x4 >= 1.) {
		if (x3 < 1.) {
			x3 = 1.
			xerr = 1
		}
		if (x4 > xmaxf) {
			x4 = xmaxf
			xerr = 1
		}
		x1 = max (x1, 1.)
		x2 = max (x2, 1.)
		x5 = min (x5, xmaxf)
		x6 = min (x6, xmaxf)

		if (extr_sky == YES) {
# Lower window
			i1 = int (x1)
			i2 = int (x2)
			if (real (i2) >= x2) {		# to avoid going over
				i2 = i2 - 1		# due to roundoff
			}
# sum up the whole pixels
			do i = i1, i2
				skyval = skyval + slit[i]
# remove partial edges
			skyval = skyval - slit[i1] * (x1 - i1)
			skyval = skyval - slit[i2] * (i2 + 1 - x2)

# Upper window
			i1 = int (x5)
			i2 = int (x6)
			do i = i1, i2
				skyval = skyval + slit[i]
			skyval = skyval - slit[i1] * (x5 - i1)
			skyval = skyval - slit[i2] * (i2 + 1 - x6)
		}

# Object window
		i1 = int (x3)
		i2 = int (x4)
		do i = i1, i2
			objval = objval + slit[i]
		objval = objval - slit[i1] * (x3 - i1)
		objval = objval - slit[i2] * (i2 + 1 - x4)

# Scale sky (if extracted):
		if (extr_sky == YES)
			skyval = skyval * (x4 - x3) / (x2 - x1 + x6 - x5)
	} else {
		xerr = 1
	}

	return (xerr)
end

######### DEBUG_MODULE

procedure	debug_mod (e1, a2, a3, a4, ccd, sys)

double	e1[3,3], a2[3,3], a3[3,3], a4[3,3]	# transforms
real	ccd[NCCD,3]				# CCD geometry
double	sys[NPARAM]				# system parameters

double	r[3]				# direction vector
double	alpha, beta, gamma		# grating angles

double	xslit, yslit, w
int	n
real	xccd, yccd
real	xics, yics
real	x1, x2, y1, y2

begin
	xccd = 1024
	yccd = 2048
	n = 2
	xccd = 1024
	yccd = 2048
	n = 4

	xslit = -74.2
	yslit = 68.026
	w = 0.696543

	call pt_xfm (xslit, yslit, w, e1, a2, a3, a4, ccd, sys, xccd, yccd, n)
call eprintf ("DEBUG_MOD: PT_XFM\n")
	call eprintf ("xy: %7.2f,%7.2f  chip=%d \n")
call pargr (xccd)
call pargr (yccd)
call pargi (n)

	call pre_grating (xslit, yslit, e1, a2, sys, r)

# transform into the grating system
	call gen_xfm (r, a3, YES)

	alpha = -atan2 (-r[2], -r[3])
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))
call eprintf ("gamma = %7.4f PRE_GRATING (OK)\n\n")
call pargd (RADTODEG(gamma))

	beta = asin ((ORDER(sys) * GRLINES(sys) * w / cos (gamma)) - sin (alpha))
	call grat_to_ics (beta, gamma, sys, a4, a3, xics, yics)
	call ics_to_ccd (xics, yics, ccd, n, x1, y1)
call eprintf ("DEBUG_MOD: replicate pt_xfm: beta=%7f  x,y= %7f,%7f,%2d\n\n")
call pargd (RADTODEG(beta))
call pargr (x1)
call pargr (y1)
call pargi (n)




	call ccd_to_ics (xccd, yccd, ccd, n, xics, yics)
#	n = ident_ccd (xpix, ypix)
	call ics_to_ccd (xics, yics, ccd, n, x1, y1)

	call eprintf ("DEBUG MODULE: ICS\n")
	call eprintf ("xy: %7.2f,%7.2f  %7.2f,%7.2f   %7.2f,%7.2f\n\n")
		call pargr (xccd)
		call pargr (yccd)
		call pargr (xics)
		call pargr (yics)
		call pargr (x1)
		call pargr (y1)


	call ics_to_grat (xics, yics, sys, a4, a3, beta, gamma)
	call grat_to_ics (beta, gamma, sys, a4, a3, x2, y2)
	
	call eprintf ("DEBUG MODULE: ICS-GRAT \n")
	call eprintf ("xy: %7.2f,%7.2f  %7.4f,%7.4f   %7.2f,%7.2f\n\n")
		call pargr (xics)
		call pargr (yics)
		call pargd (RADTODEG(beta))
		call pargd (RADTODEG(gamma))
		call pargr (x2)
		call pargr (y2)

call eprintf ("Apparent large round-off (>1.e-4) in x above ... \n\n")

		call ccd_to_grating (xccd, yccd, n, sys, ccd, a4, r)

		call gen_xfm (r, a3, YES)

		beta = -atan2 (r[2], r[3])
		gamma = asin (r[1])

	call eprintf ("DEBUG MODULE: CCD_GRAT \n")
	call eprintf ("xy: %7.2f,%7.2f  %7.4f,%7.4f  \n\n")
		call pargr (xics)
		call pargr (yics)
		call pargd (RADTODEG(beta))
		call pargd (RADTODEG(gamma))


call eprintf ("DEBUG gen_xfm, r= %12.9f %12.9f %12.9f err=%7g\n")
call pargd (r[1])
call pargd (r[2])
call pargd (r[3])
call pargd (r[1]*r[1]+r[2]*r[2]+r[3]*r[3] - 1.0D0)
		call gen_xfm (r, a3, NO)
		call gen_xfm (r, a3, YES)
call eprintf ("DEBUG gen_xfm, r= %12.9f %12.9f %12.9f err=%7g\n")
call pargd (r[1])
call pargd (r[2])
call pargd (r[3])
call pargd (r[1]*r[1]+r[2]*r[2]+r[3]*r[3] - 1.0D0)

# call fatal (0, "TEST_MODULE terminate")
end

######### TEST_MODULE: Efficiency

procedure	blaze_mod (alpha0, sys, lamb, lamr)

double	alpha0					# input
double	sys[NPARAM]				# system parameters
double	lamb, lamr

double	beta0			# input beta
double	gamma0, cosg		# input gamma

double	wave			# wavelength
double	alpha, beta		# ghost grating angles
double	arg

double	gamma		# gamma param

double	delta, cosd
real	sum
real	bf[5]

int	i, j

begin
	delta = DEGTORAD(26.1)
	cosd = cos (delta)
	gamma0 = 0.
	cosg = cos (gamma0)

	do j = 0, 40 {
		wave = lamb + (lamr-lamb) * (j / 40.)
		beta0 = asin ((ORDER(sys) * GRLINES(sys) * wave / cosg) - sin (alpha0))

		alpha = beta0
		sum = 0.
		do i = -2, 2 {
			arg = (i * GRLINES(sys) * wave / cosg) - sin (alpha)
			if (arg < -1. || arg > 1.) {
				bf[i+2+1] = 0.
				next
			}

			beta = asin (arg)

			gamma = PI / GRLINES(sys) * cosd / wave *
				(sin (alpha - delta) + sin (beta - delta))

			if (abs (gamma) < 1.d-9) {
			    bf[i+2+1] = 1.
			} else {
			    bf[i+2+1] = sin (gamma) ** 2 / (gamma * gamma)
			    sum = sum + bf[i+2+1]
			}
		}
		call printf ("%5f %5f   %9.7f %9.7f %9.7f %9.7f %9.7f\n")
			call pargd (wave)
			call pargr (1.-sum)
			call pargr (bf[1])
			call pargr (bf[2])
			call pargr (bf[3])
			call pargr (bf[4])
			call pargr (bf[5])
	}

call flush (STDOUT)
call fatal (0, "BLAZE_MODULE terminate")
end
