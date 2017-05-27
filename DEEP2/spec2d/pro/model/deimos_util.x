# NOTE that proj_to_mask, mask_to_proj should be rewritten to allow inputs
# same as outputs!

# DEIMOS_UTIL:  various utilities for DEIMOS,
#	PROJ_TO_MASK: project planar onto curved slitmask coordinates
#	MASK_TO_PROJ: project slimask coords (curved surface) onto plane
# XXX	CCD_MAP:      Map slitmask onto CCD coords (revised)
# XXX	SKY_GEOM:     get various angles on sky (LRIS-specific in defn of PA)

include	<math.h>
include <math/gsurfit.h>
include	"deimos.h"
include	"keck.h"

#
# PROJ_TO_MASK: project planar onto curved slitmask coordinates
# Double inputs, outputs
# Note that this is pure geometry -- any empirically determined corrections
# should go elsewhere...
#

procedure proj_to_mask (xp, yp, ap, xc, yc, ac)

double	xp, yp			# x,y in (focal)-plane system
double	ap			# position angle in planar system (deg)
double	xc, yc			# returned x,y values on mask surface
double	ac			# returned angle on curved surface (deg)

double	mu, cosm			# mu, cos(mu)
double	cost, tant			# cos, tan of mask tilt angle
double	tanpa				# tan PA

double	rho			# radius from telescope optical axis
double	hs, hm			# height of image surface, mask above datum
double	xx, yy		# Work variables corresponding to xc, yc

begin
	mu = asin (xp / M_RCURV)
	cosm = cos (mu)
	cost = cos (DEGTORAD(M_ANGLE))
	tant = tan (DEGTORAD(M_ANGLE))
	xx =  M_RCURV * mu
	yy =  (yp - ZPT_YM) / cost + M_RCURV * tant * (1. - cosm)

	tanpa = tan (DEGTORAD(ap)) * cosm / cost + tant * xp / M_RCURV
	ac = RADTODEG(atan (tanpa))

# What follows is a small correction for the fact that the mask does
# not lie exactly in the spherical image surface (where the distortion
# values and gnomonic projection are calculated) and the rays are arriving
# from the pupil image; thus, the exact locations are moved _slightly_
# wrt the telescope optical axis.  Note also that these corrections are
# only calculated to first order.

# Spherical image surface height:
	rho = sqrt (xp * xp + yp * yp)
	hs = R_IMSURF * (1. - sqrt (1. - (rho / R_IMSURF) ** 2))
# Mask surface height:
	hm = MASK_HT0 + yy * sin (DEGTORAD(M_ANGLE)) + M_RCURV * (1. - cosm)
# Correction:
	yc = yy + (hs - hm) * yp / PPLDIST / cost
	xc = xx + (hs - hm) * xp / PPLDIST / cosm
end

#
# MASK_TO_PROJ: project slitmask coords (curved surface) onto plane
# Double inputs, outputs
# Note that this is pure geometry -- any empirically determined corrections
# should go elsewhere...
#

procedure mask_to_proj (xc, yc, ac, xp, yp, ap)

double	xc, yc			# x,y values on mask surface
double	ac			# position angle on curved surface
double	xp, yp			# returned x,y in (focal)-plane system
double	ap			# returned position angle in planar system

double	mu, cosm			# mu, cos (mu)
double	cost, tant			# cos, tan of mask tilt angle
double	tanpa				# tan PA

double	rho			# radius from telescope optical axis
double	hs, hm			# height of image surface, mask above datum
double	xx, yy		# Work variables corresponding to xp, yp

begin
	mu = xc / M_RCURV
	cosm = cos (mu)
	cost = cos (DEGTORAD(M_ANGLE))
	tant = tan (DEGTORAD(M_ANGLE))
	xx =  M_RCURV * sin (mu)
	yy =  (yc - M_RCURV * tant * (1. - cosm)) * cost + ZPT_YM

	tanpa = (tan (DEGTORAD(ac)) - tant * xx / M_RCURV) * cost / cosm
	ap = atan (tanpa)

# What follows is a small correction for the fact that the mask does
# not lie exactly in the spherical image surface (where the distortion
# values and gnomonic projection are calculated) and the rays are arriving
# from the pupil image; thus, the exact locations are moved _slightly_
# wrt the telescope optical axis.  Note also that these corrections are
# only calculated to first order.

# Spherical image surface height:
	rho = sqrt (xx * xx + yy * yy)
	hs = R_IMSURF * (1. - sqrt (1. - (rho / R_IMSURF) ** 2))
# Mask surface height:
	hm = MASK_HT0 + yc * sin (DEGTORAD(M_ANGLE)) + M_RCURV * (1. - cosm)
	yp = yy - (hs - hm) * yy / PPLDIST 
	xp = xx - (hs - hm) * xx / PPLDIST 
end



#########################################################
#
# CAM_DISTORT: remove camera distortion from distorted angles.
## Currently DEIMOS theoretical distortion curves apply. !!!
# Double inputs, outputs; bool apply
#
 
# These coeff's from polyfit; for legendre, must translate the coeffs:
#define		CAMD_C0		 0.0013801
#define		CAMD_C2		 0.0578925
#define		CAMD_C4		-1.144587

#define		CAMD_C0		 0.			# should be 0!
#define		CAMD_C2		 0.5912466		# adjust for scale
#define		CAMD_C4		 0.

define		MAX_TAN		 0.07852	# tan (angle) at one side

#     1.001379   1.374380E-9    0.05754194  -1.423267E-7      -1.18797
#      1.       -2.52929E-10     0.0457563   2.308757E-8    -0.3088123
#    -4.616165E-7       -14.917
define		CAMD_C0		 1.D0
define		CAMD_C2		 0.0457563
define		CAMD_C4		 -0.3088123
define		CAMD_C6		 -14.917

procedure	cam_distort (thetain, thetaout, apply)

double	thetain, thetaout		# angles in, out
bool	apply				# apply (vs. remove) distortion?

begin
	if (apply) {
		thetaout = thetain *
			(CAMD_C0 + thetain*thetain *
			(CAMD_C2 + thetain*thetain *
			(CAMD_C4 + thetain*thetain * CAMD_C6)))
	} else {
		thetaout = thetain /
			(CAMD_C0 + thetain*thetain *
			(CAMD_C2 + thetain*thetain *
			(CAMD_C4 + thetain*thetain * CAMD_C6)))
	}

end

#
# CCD_MAP: map x,y (slitmask) into CCD coords. Assumes fit from geomap entered
# in def_map()

# To produce the coeffs, run geomap on a file that has X,Y,x,y (SM,CCD).
# The two columns of coeffs are for the x-fit and y-fit respectively.
# Note that the linear terms (surface1) must be added to the first 3
# coeffs of the distortion terms (surface2)
#
# Mod 15-mar-99: add call to inst_config to get map name.

procedure	ccd_map (xmm, ymm, xccd, yccd, npts)

real	xmm[npts], ymm[npts]			# input slitmask coords
real	xccd[npts], yccd[npts]			# output ccd coords
int	npts					# number of slits

char	input[SZ_FNAME]				# name of map
pointer	fd					# file descriptor
pointer	sfx, sfy

pointer	open()
begin
#	call gs_ingest ("newlris$mappings/lris.sm2ccd", sfx, sfy)
	call clgstr ("inst_config.sm2ccd", input, SZ_FNAME)
	fd = open (input, READ_ONLY, TEXT_FILE)
	call gs_ingest (fd, sfx, sfy)
	call gsvector (sfx, xmm, ymm, xccd, npts)
	call gsvector (sfy, xmm, ymm, yccd, npts)
	call gsfree (sfx)
	call gsfree (sfy)
	call close (fd)

end

# DEF_CCD_MAP: define the slit-ccd mapping  OBSOLETE?

procedure	def_ccd_map (sfx, sfy)

# These coefficients were taken from a fit to the M71_D1 map.
# Additional offsets correct to LRIS stow position (zeropt for flexure mapping)

pointer	sfx, sfy			# pointers to surface fits in x,y

int	ncoeff
pointer	xcoeff, ycoeff			# coeff's in x,y

begin
	ncoeff = 20
	call malloc (xcoeff, ncoeff, TY_REAL)
	call malloc (ycoeff, ncoeff, TY_REAL)

	Memr[xcoeff  ]  = 2.
	Memr[xcoeff+1]  = 4.
	Memr[xcoeff+2]  = 3.
	Memr[xcoeff+3]  = 1.
	Memr[xcoeff+4]  = 1.
	Memr[xcoeff+5]  = 216.
	Memr[xcoeff+6]  = 1.
	Memr[xcoeff+7]  = 336.
	Memr[xcoeff+8]  =  0.4259691282581062 + 1089.395693397732 - 42. - 0.7 - 18.3		# -14.4 for mirror realignment
	Memr[xcoeff+9]  =  1.135406932277368  + 694.8209734223686
	Memr[xcoeff+10] =  0.4332245141494603
	Memr[xcoeff+11] =  0.6601396941753989
	Memr[xcoeff+12] = -0.364762294701932 + 2.9168802655806751
	Memr[xcoeff+13] = -1.114474752798539
	Memr[xcoeff+14] = -0.07100544020606963
	Memr[xcoeff+15] =  0.01848460943466597
	Memr[xcoeff+16] =  1.418819558933747
	Memr[xcoeff+17] =  3.289284015335439
	Memr[xcoeff+18] =  0.08722518142092386
	Memr[xcoeff+19] = -0.06158022096394576

	Memr[ycoeff  ]  = 2.
	Memr[ycoeff+1]  = 3.
	Memr[ycoeff+2]  = 4.
	Memr[ycoeff+3]  = 1.
	Memr[ycoeff+4]  = 1.
	Memr[ycoeff+5]  = 216.
	Memr[ycoeff+6]  = 1.
	Memr[ycoeff+7]  = 336.
	Memr[ycoeff+8]  =  0.464683175237381  + 1096.128797371226  - 6.0 - 70.3			# -67.2 for mirror realignment
	Memr[ycoeff+9]  =  0.1765672348323415 +  1.460280443638984
	Memr[ycoeff+10] =  0.1613619533305995
	Memr[ycoeff+11] = -2.080077948000524 + -1093.182286712313
	Memr[ycoeff+12] = -3.252201997220546
	Memr[ycoeff+13] = -1.548886149883428
	Memr[ycoeff+14] =  1.719400498669824
	Memr[ycoeff+15] = -0.02326077701339106
	Memr[ycoeff+16] =  0.3049788335944159
	Memr[ycoeff+17] = -2.430251103657301
	Memr[ycoeff+18] = -0.4204804794749471
	Memr[ycoeff+19] =  0.2496466944483961

	call gsrestore (sfx, Memr[xcoeff], ncoeff)
	call gsrestore (sfy, Memr[ycoeff], ncoeff)

	call mfree (xcoeff, TY_REAL)
	call mfree (ycoeff, TY_REAL)
end

#
# GS_INGEST: Read in a geomap surface fit file and add the two surfaces.
# Note: This procedure may be easily modified to apply a rotation/scale/offset
# via modification of the FIRST surface. Alternatively, the input file may be
# modified.
#

procedure	gs_ingest (fd, gsx, gsy)

pointer	gsx, gsy				# The returned gsurfit ptrs
pointer	fd					# file descriptor

char	label[SZ_LINE]
int	nval
double	arg1, arg2
int	i
pointer	sp
pointer	bufx, bufy
pointer	gsx1, gsy1, gsx2, gsy2		# pointers to separated surfaces

bool	streq()
int	fscan(), nscan()

begin
# call eprintf ("Open map file: %s\n")
# call pargstr (input)
#	if (access (input, READ_ONLY, TEXT_FILE) == YES) {
#		fd = open (input, READ_ONLY, TEXT_FILE)
#	} else {
#		call fatal (0, "Cannot open surface mapping file! \n")
#	}

	call seek (fd, BOF)

	while (fscan (fd) != EOF) {
		call gargwd (label, SZ_LINE)
		if (streq (label, "surface1")) {
			call gargi (nval)
			call smark (sp)
			call salloc (bufx, nval, TY_REAL)
			call salloc (bufy, nval, TY_REAL)

			do i = 0, nval-1 {
				if (fscan (fd) == EOF)
				    call fatal (0, "Unexpected EOF in map!\n")
				call gargd (arg1)
				call gargd (arg2)
				if (nscan() < 2)
				    call fatal (0, "Incorrect surface map!\n")
				Memr[bufx+i] = arg1
				Memr[bufy+i] = arg2
			}
# restore the fit (don't know if nval is needed...)
			call gsrestore (gsx1, Memr[bufx], nval)
			call gsrestore (gsy1, Memr[bufy], nval)
			call sfree (sp)
		}
		if (streq (label, "surface2")) {
			call gargi (nval)
			call smark (sp)
			call salloc (bufx, nval, TY_REAL)
			call salloc (bufy, nval, TY_REAL)

			do i = 0, nval-1 {
				if (fscan (fd) == EOF)
				    call fatal (0, "Unexpected EOF in map!\n")
				call gargd (arg1)
				call gargd (arg2)
				if (nscan() < 2)
				    call fatal (0, "Incorrect surface map!\n")
				Memr[bufx+i] = arg1
				Memr[bufy+i] = arg2
			}
# restore the fit (don't know if nval is needed...)
			call gsrestore (gsx2, Memr[bufx], nval)
			call gsrestore (gsy2, Memr[bufy], nval)
			call sfree (sp)
		}
	}

#	call close (fd)

# Add the two mappings
	call gsadd (gsx1, gsx2, gsx)
	call gsadd (gsy1, gsy2, gsy)

	call gsfree (gsx1)
	call gsfree (gsy1)
	call gsfree (gsx2)
	call gsfree (gsy2)
end

#
# ZETA_MAP: map x,y (slitmask) into ZETA angles. Assumes gs_ingest already run.
# To produce the coeffs, run grid_map in "zeta_mode", then run geomap on result.
#

procedure	zeta_map (sfx, sfy, xmm, ymm, zetax, zetay, npts)

pointer	sfx, sfy				# pointers to fits
real	xmm[npts], ymm[npts]			# input slitmask coords
real	zetax[npts], zetay[npts]		# output ccd coords
int	npts					# number of slits

begin
	call gsvector (sfx, xmm, ymm, zetax, npts)
	call gsvector (sfy, xmm, ymm, zetay, npts)
end



#	
# SKY_GEOM: get mask, dispersion PAs on sky
#

procedure	sky_geom (im, maskpa, z, adpa)

pointer	im				# image pointer
real	maskpa				# mask PA (deg)
real	z				# zenith distance (deg)
real	adpa				# atm. disp. PA (deg) blue to red

real	w1
real	dec, ha, az, atmdisp, expo

int	imaccf()
real	imgetr()
begin
	if (imaccf (im, "ROTPOSN") == YES) {
		maskpa = imgetr (im, "ROTPOSN") + 90.
	} else {
		call eprintf ("No mask angle found -- assume 0 \n")
		maskpa = 0.
	}

	if (imaccf (im, "DEC") == YES && imaccf (im, "HA") == YES) {
		dec  = imgetr (im, "DEC")
		dec = min (dec, 89.5)
		ha   = imgetr (im, "HA")
		expo = imgetr (im, "EXPOSURE")
		ha = ha + 0.5*expo/3600.
		w1 = 6000.			# dummy here
		call atm_geom (ha, dec, w1, w1, OBS_LAT, z, az, adpa, atmdisp)
call eprintf ("SKY_GEOM:  maskpa, adpa, z (az): %5.1f %5.1f %4.1f (%5.1f)\n")
call pargr (maskpa)
call pargr (adpa)
call pargr (z)
call pargr (az)
	} else {
		call eprintf ("WARNING: No Elevation -- 90 assumed\n")
		z = 0.
		adpa = 0.
	}
end
