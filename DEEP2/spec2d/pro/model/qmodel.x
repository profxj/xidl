# QMODEL: simplified model; calls in the prev. gen maps and grating pars

include	<math.h>
include	<math/gsurfit.h>
include	"deimos.h"
include	"instrument.h"


procedure	t_qmodel()

char	amap[SZ_FNAME], bmap[SZ_FNAME]		# input mappings
real	xmm, ymm, wave			# X,Y in slitmask, wave
pointer	fda, fdb

double	sys[NPARAM]				# system parameters
double	a3[3,3]				# grating transform
double	r[3]
double	alpha, beta, gamma
real	tanx, tany			# should be double; check mapping evals
real	xics, yics			# pixel values in ICS
real	xpix, ypix			# pixel values on CCD
real	cost, sint

pointer	asfx, asfy			# pointers to surface fits (amap)
pointer	bsfx, bsfy			# pointers to surface fits (bmap)

int	clgeti()
real	clgetr()
real	gseval()
pointer	open()

begin
	xmm = clgetr ("xmm")
	ymm = clgetr ("ymm")
	wave = 1.e-4 * clgetr ("wave")			# in microns

	MU(sys) = DEGTORAD (clgetr ("mu"))
	GR_YERR(sys) = DEGTORAD (clgetr ("roll3"))
	GR_ZERR(sys) = DEGTORAD (clgetr ("o3"))

	ORDER(sys) = clgeti ("norder")
	GRLINES(sys) = 1.e-3 * clgeti ("gmm")		# in microns

	call clgstr ("amap", amap, SZ_FNAME)
	fda = open (amap, READ_ONLY, TEXT_FILE)

	call clgstr ("bmap", bmap, SZ_FNAME)
	fdb = open (bmap, READ_ONLY, TEXT_FILE)

# Iitialize the maps
	call gs_ingest (fda, asfx, asfy)
	call gs_ingest (fdb, bsfx, bsfy)

#######################################################################
# Below is the actual calculation; gseval simply evaluates the mappings
#######################################################################

# set up the grating transform
	call gsetup (a3, sys)

# Get mapping and convert to r[3]
	tanx = gseval (asfx, xmm, ymm)
	tany = gseval (asfy, xmm, ymm)

	r[3] = -1.d0 / sqrt (1. + tanx*tanx + tany*tany)
	r[1] = r[3] * tanx
	r[2] = r[3] * tany

# xform into grating system
	call gen_xfm (r, a3, YES)

# convert to alpha,gamma
	alpha = -atan2 (-r[2], -r[3])
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))

# Aplly the grating equation
	beta = asin ((ORDER(sys)*GRLINES(sys)*wave / cos (gamma)) - sin (alpha))

# convert beta, gamma into x,y,z (cf Schroeder p259); note sign reversal of beta
	r[1] = sin (gamma)
	r[2] = sin (-beta) * cos (gamma)
	r[3] = cos (-beta) * cos (gamma)

# xform out of grating system
	call gen_xfm (r, a3, NO)

# convert to tanx, tany
	tanx = (-r[1] / -r[3])
	tany = (-r[2] / -r[3])

# get mapping into ICS pixels
	xics = gseval (bsfx, tanx, tany)
	yics = gseval (bsfy, tanx, tany)

# get mapping into CCD pixels
	cost = cos (DEGTORAD(90.-0.27))
	sint = sin (DEGTORAD(90.-0.27))
# the offset should be: 
#	2. * CCDXPIX + 4.*CCDXEDG/PIX_SZ + 1.5*NOMXGAP/PIX_SZ +
#		FCSYEDG/PIX_SZ * 0.5*FCSYPIX + gap = 4546.4 + gap
	xics = xics - 4584.		# (relative to center of FCS2)
	yics = yics - -600. 		# (relative to center of FCS2)
	xpix =  xics * cost - yics * sint + 600.
	ypix = -xics * sint - yics * cost + 300.

call eprintf ("%8.3f %7.3f %8.2f -->ICS: %7.1f %7.1f -->FC2: %6.1f %6.1f \n")
	call pargr (xmm)
	call pargr (ymm)
	call pargr (wave*1.e4)
	call pargr (xics)
	call pargr (yics)
	call pargr (xpix)
	call pargr (ypix)

end



procedure	gsetup (a3, sys)

double	a3[3,3]			# grating transforms
double	sys[NPARAM]		# system parameters

double	thetan
double	cost, sint
double	xsi, cosx, sinx
double	rhon
double	cosr, sinr

begin
# ... below assumes phin=0. (ie adopts the roll/yaw approach)

	thetan = -MU(sys)
	xsi = GR_ZERR(sys)
	rhon = GR_YERR(sys)
	cost = cos (thetan)
	sint = sin (thetan)
	cosx = cos (xsi)
	sinx = sin (xsi)
	cosr = cos (rhon)
	sinr = sin (rhon)

	a3[1,1] =  cosx*cosr
	a3[1,2] =  sint*sinr + cost*sinx*cosr
	a3[1,3] = -cost*sinr + sint*sinx*cosr
	a3[2,1] = -sinx
	a3[2,2] =  cost*cosx
	a3[2,3] =  sint*cosx
	a3[3,1] =  cosx*sinr
	a3[3,2] = -sint*cosr + cost*sinx*sinr
	a3[3,3] =  cost*cosr + sint*sinx*sinr
end



##
## GEN_XFM: general transform of r[3] into another CS desribed by a; "forward"
## is YES/NO to describe if xform is into or ou-of CS
## Note that the appropriate operation (eg transmission, reflection) must be
## applied afterward
##

#procedure	gen_xfm (r, a, forward)

#double	r[3]
#double	a[3,3]
#int	forward

#double	rp[3]
#int	i, j

#begin
## transform
#	if (forward == YES) {
#	    do i = 1, 3 {
#		rp[i] = 0.
#		do j = 1, 3 {
#			rp[i] = rp[i] + a[i,j] * r[j]
#		}
#	    }
#
#	    do i = 1, 3
#		r[i] = rp[i]
#
#	} else {
#		
#	    do i = 1, 3
#		rp[i] = r[i]
#
#	    do i = 1, 3 {
#		r[i] = 0.
#		do j = 1, 3 {
#			r[i] = r[i] + a[j,i] * rp[j]
#		}
#	    }
#	}
#end
