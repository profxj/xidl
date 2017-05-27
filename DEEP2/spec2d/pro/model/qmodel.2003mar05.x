
# QMODEL: simplified model; calls in the prev. gen maps and grating pars

include	<math.h>
include	<math/gsurfit.h>
include	"deimos.h"
include	"instrument.h"


procedure	t_qmodel()

char	amap[SZ_FNAME], bmap[SZ_FNAME]		# input mappings
char	input[SZ_FNAME]				# input file of x,y pairs
real	xmm, ymm, wave			# X,Y in slitmask, wave
real	scaling				# scaling difference from 1
pointer	fd

pointer	fda, fdb
pointer	map[4]				# pointers to surf fits (1,2=amap;3,4=b)

double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]			# CCD geometry
double	a3[3,3]				# grating transform
real	xics, yics			# pixel values in ICS
real	xpix, ypix			# pixel values on CCD
int	stat
int	chip

bool	file_list			# is there a file with x,y list?

char	id[16]		# TMP

int	qxfm(), qxfm_init()

bool	strne()
int	clgeti()
int	fscan(), nscan()
real	clgetr()
pointer	open()

begin
call eprintf ("INDV MODE PARTIALLT DISABLED!\n")

	call clgstr ("input", input, SZ_FNAME)
	file_list = strne (input, "")
	if (file_list) {
		fd = open (input, READ_ONLY, TEXT_FILE)
	} else {
		xmm = clgetr ("xmm")
		ymm = clgetr ("ymm")
	}

	wave = 1.e-4 * clgetr ("wave")			# in microns
	scaling = 1. + clgetr ("scale_adj")

	MU(sys) = DEGTORAD (clgetr ("mu"))
	GR_YERR(sys) = DEGTORAD (clgetr ("roll3"))
	GR_ZERR(sys) = DEGTORAD (clgetr ("o3"))

	ORDER(sys) = clgeti ("norder")
	GRLINES(sys) = 1.e-3 * clgeti ("gmm")		# in microns

	call clgstr ("amap", amap, SZ_FNAME)
	fda = open (amap, READ_ONLY, TEXT_FILE)

	call clgstr ("bmap", bmap, SZ_FNAME)
	fdb = open (bmap, READ_ONLY, TEXT_FILE)

# Initialize the mappings
	stat = qxfm_init (fda, fdb, map, a3, sys, ccd)

# Loop if need be
	while (fscan (fd) != EOF) {
		call gargr (xmm)
		call gargr (ymm)
		call gargstr (id,16)	# TMP
		if (nscan() < 2)
			next

# calculate the mapping
	stat = qxfm (map, a3, sys, ccd, xmm, ymm, wave, scaling, xics, yics, xpix, ypix, chip, YES, YES)

	if (stat != ON_CHIP) {
call eprintf ("%8.3f %7.3f %8.2f -->ICS: %7.1f %7.1f -->OFF CHIPS \n")
		call pargr (xmm)
		call pargr (ymm)
		call pargr (wave*1.e4)
		call pargr (xics)
		call pargr (yics)
call printf ("%s OFF_CHIP\n")
call pargstr (id)		# TMP
	} else {
call eprintf ("%8.3f %7.3f %8.2f -->ICS: %7.1f %7.1f --> %6.1f %6.1f (%d)\n")
		call pargr (xmm)
		call pargr (ymm)
		call pargr (wave*1.e4)
		call pargr (xics)
		call pargr (yics)
		call pargr (xpix)
		call pargr (ypix)
		call pargi (chip)
call printf ("%6.1f %6.1f %-16s\n")
call pargr (xpix)
call pargr (ypix)
call pargstr (id)		# TMP
	}
	}

	if (file_list)
		call close (fd)

end




int	procedure qxfm (map, a3, sys, ccd, xmm, ymm, wave, scaling, xics, yics, xpix, ypix, chip, find_chip, mos_coords)

pointer	fda, fdb		# pointers to map input files

pointer	map[4]			# pointer to maps
double	a3[3,3]			# grating transform
double	sys[NPARAM]		# system parameters
real	ccd[NCCD,3]		# CCD geometry
real	xmm, ymm		# x,y slitmask coords
real	wave			# wavelength (um)
real	scaling			# scaling adjustment
real	xics, yics		# pixel values in ICS
real	xpix, ypix		# pixel values on CCD (PANE coordinates)
int	chip			# CCD number
int	stat			# on/off chip status
int	find_chip		# find the chip (or force to existing value)?
int	mos_coords		# return mosaic (or individual) coords?

double	r[3]
double	alpha, beta, gamma
real	tanx, tany		# should be double; check mapping evals
int	n

int	ident_ccd()
real	gseval()

int	qxfm_init()
begin

# Get mapping and convert to r[3]
	tanx = gseval (map[1], xmm, ymm)
	tany = gseval (map[2], xmm, ymm)
# call eprintf ("tanx,y: %5f %5f\n")
# call pargr (tanx)
# call pargr (tany)

	r[3] = -1.d0 / sqrt (1. + tanx*tanx + tany*tany)
	r[1] = r[3] * tanx
	r[2] = r[3] * tany

# xform into grating system
	call gen_xfm (r, a3, YES)

# convert to alpha,gamma
	alpha = -atan2 (-r[2], -r[3])
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))

# Apply the grating equation
	beta = asin ((ORDER(sys)*GRLINES(sys)*abs(wave) / cos (gamma)) - sin (alpha))

# convert beta, gamma into x,y,z (cf Schroeder p259); note sign reversal of beta
	if (wave > 0) {
		r[1] = sin (gamma)
		r[2] = sin (-beta) * cos (gamma)
		r[3] = cos (-beta) * cos (gamma)
	} else {
		r[1] = sin (gamma)
		r[2] = sin (beta) * cos (gamma)
		r[3] = cos (beta) * cos (gamma)
	}

# xform out of grating system
	call gen_xfm (r, a3, NO)

# convert to tanx, tany
	tanx = (-r[1] / -r[3])
	tany = (-r[2] / -r[3])

# get mapping into ICS pixels
	xics = gseval (map[3], tanx, tany)
	yics = gseval (map[4], tanx, tany)

	xics = xics * scaling
	yics = yics * scaling

	n = ident_ccd (xics, yics)

	if (find_chip == YES) {
		chip = n
	} else {
		n = chip		# XXX check?
	}

# Put coordinates into desired form:
	if (mos_coords == YES) {
		call ics_to_ccd (xics, yics, ccd, n, xpix, ypix, stat)
# Convert to full mosaic image
		call mosim_coord (xpix, ypix, n)
	} else {
		call ics_to_ccd (xics, yics, ccd, chip, xpix, ypix, stat)
	}
	

#call eprintf ("%8.3f %7.3f %8.2f -->ICS: %7.1f %7.1f -->FC2: %6.1f %6.1f \n")
	return (stat)


entry	 qxfm_init (fda, fdb, map, a3, sys, ccd)

# Iitialize the maps
	call gs_ingest (fda, map[1], map[2])
	call gs_ingest (fdb, map[3], map[4])

# set up the grating transform
	call gsetup (a3, sys)

# CCD: define the geometry of the mosaic
	call ccd_geom (ccd, sys, YES)
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
