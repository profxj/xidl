# TR_UTIL: package of tracing routines
include	<math.h>

procedure	reflect (inx, iny, tnx, tny, erx, ery, rfx, rfy)

# Note that eventually, all should be doubles

double	inx, iny			# input angles
double	tnx, tny			# angle of normals
double	erx, ery			# alignment errors
double	rfx, rfy			# output angles

double	x, y, z, xp, yp, zp
double	tanx, tany
double	theta, phi
double	cost, sint, cosp, sinp

begin

# rotate to coord system of normal
	tanx = tan (tnx + erx)
	tany = tan (tny + ery)

	phi = atan2 (tany, tanx) + HALFPI
	theta = atan (sqrt (tanx*tanx + tany*tany))
	cosp = cos (phi)
	sinp = sin (phi)
	cost = cos (theta)
	sint = sin (theta)

# allow psi to be zero

# Now consider the (x,y,z) of the incoming ray; take z to be unity
	z = 1
	x = z * tan (inx)
	y = z * tan (iny)

# transform to primed system
	xp =       cosp * x  +       sinp * y
	yp = -cost*sinp * x  +  cost*cosp * y  +  sint * z
	zp =  sint*sinp * x  -  sint*cosp * y  +  cost * z

# upon reflection,
	xp =  xp
	yp =  yp
	zp = -zp	# have potentially mixed some sign conventions here

# and the reverse transform:
	x = cosp * xp  -  cost*sinp * yp  +  sint*sinp * zp
	y = sinp * xp  +  cost*cosp * yp  -  sint*cosp * zp
	z =                    sint * yp  +       cost * zp

# and
	rfx = atan (x/z)
	rfy = atan (y/z)
end





# PBOLA_XY:  Get the height and coordinates at a parabolic surface.
# x0, y0 are x, y in simple projection to a plane at location of coll.
# Implicitly assumes everything is on-axis!

procedure	pbola_xy (x0, y0, tanx, tany, rcurv, xf, yf, h)

double	x0, y0		# coordinates of ray projected onto plane
double	tanx, tany	# tan of projected angles in x,y
double	rcurv		# radius of curvature of paraboloid
double	xf, yf		# (x,y) at intercept
double	h		# height at intercept

double	tant, q
double	fact		# fractional change, delr/r0

begin
	tant = sqrt (tanx*tanx + tany*tany)
	
	q = sqrt (x0*x0 + y0*y0) / rcurv * tant

	fact = 1. / q * (1. + q - sqrt (1. + 2.*q))	# delr / r0
	xf = x0 * (1.- fact)
	yf = y0 * (1.- fact)
	h = (xf*xf + yf*yf) / (2 * rcurv)
end





# SPHER_XY:  Get the height and coordinates at a spherical surface.
# x0, y0 are x, y in simple projection to a plane at location of coll.
# Implicitly assumes everything is on-axis!

procedure	spher_xy (x0, y0, tanx, tany, rcurv, xf, yf, h)

double	x0, y0		# coordinates of ray projected onto plane
double	tanx, tany	# tan of projected angles in x,y
double	rcurv		# radius of curvature of sphere
double	xf, yf		# (x,y) at intercept
double	h		# height at intercept

double	tant, r0

begin
	tant = sqrt (tanx*tanx + tany*tany)
	r0 = sqrt (x0*x0 + y0*y0)

# solve for height
	h = (r0*tant + rcurv - sqrt (rcurv*rcurv + 2*rcurv*r0*tant - r0*r0)) /
							(1. + tant*tant)

	xf = x0 - h * tanx
	yf = y0 - h * tany
end


define	R_COLL	4000.D0
# define	PPLDIST	19947.D0
define	PPL_FOC 	20023.15D0
define	PPL_COLL	22220.35D0

procedure	lris_coll (x, y, ax, ay, erx, ery, rx, ry, h)

real	x, y
double	ax, ay
double	erx, ery
double	rx, ry, h

double	tnx, tny

begin
	tnx = atan (x / R_COLL) 
	tny = atan (y / R_COLL)
	h = (x*x + y*y) / (2*R_COLL)
	
	call reflect (ax, ay, tnx, tny, erx, ery, rx, ry)
end

procedure	collimator (ax, ay, erx, ery, rx, ry, x, y, h)

double	ax, ay			# input angles in x, y-planes
double	erx, ery		# alignment errors in x, y-planes
double	rx, ry			# reflected angle
double	x, y, h			# x, y, height at point of reflection

double	tanx, tany		# tangents of input angles
double	x0, y0			# x, y at projection onto plane at collimator
double	tnx, tny		# normal angles

begin
# project angles to plane at collimator
	tanx = tan (ax)
	tany = tan (ay)
	x0 = PPL_COLL * tanx
	y0 = PPL_COLL * tany

# get actual reflection point
#	call spher_xy (x0, y0, tanx, tany, R_COLL, x, y, h)
	call pbola_xy (x0, y0, tanx, tany, R_COLL, x, y, h)

# angles of normal, height at that point:
	tnx = atan (x / R_COLL) 
	tny = atan (y / R_COLL)
	
	call reflect (ax, ay, tnx, tny, erx, ery, rx, ry)
end

procedure	t_testp


double	x, y, ax, ay
double	tnx, tny
double	ox, oy
double	h, dzero
real	clgetr()

begin
	x = clgetr ("xx")
	y = clgetr ("yy")
	ax = clgetr ("ax")
	ay = clgetr ("ay")

	dzero = 0.

	tnx = atan (x / R_COLL) 
	tny = atan (y / R_COLL)
	h = (x*x + y*y) / (2*R_COLL)
	
	call reflect (ax, ay, tnx, tny, dzero, dzero, ox, oy)

	call eprintf ("%9.6f %9.6f %9.6f   %9.6f %9.6f %9.6f\n")
		call pargd (ax)
		call pargd (ox)
		call pargd (tnx)
		call pargd (ay)
		call pargd (oy)
		call pargd (tny)
	call printf ("Foc: %11f %11f mm\n")
		call pargd (x/tan(ox) + h)
		call pargd (y/tan(oy) + h)

# Now suppose things work orthogonally:
		ox = 2 * tnx - ax
		oy = 2 * tny - ay
	call printf ("old: %11f %11f mm\n")
		call pargd (x/tan(ox) + h)
		call pargd (y/tan(oy) + h)
# We find the rays cross the axis at distinctly different distances.

#	call collimator (ax, ay, dzero, dzero, ox, oy, x, y, h)

#	call eprintf ("%9.6f %9.6f   %9.6f %9.6f      %7.2f %7.2f %6.3f\n")
#		call pargd (ax)
#		call pargd (ay)
#		call pargd (ox)
#		call pargd (oy)
#		call pargd (x)
#		call pargd (y)
#		call pargd (h)
end

