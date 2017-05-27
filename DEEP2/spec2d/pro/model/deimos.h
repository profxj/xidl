# DEIMOS.H -- specific parameters for DEIMOS  XXX NEEDS TO BE REVIEWED!!!
#
define	M_ANGLE	6.00D0		# Mask angle (tilt) in degrees
define	FL_TEL	150100.4D0	# Foc Len of KII WITH 3" PULLBACK (REF)
define	ZPT_YM	128.803		# Dist to tel.axis, in SMCS-XXX (mm) 5.071in
define	M_RCURV	2120.9D0	# Mask radius of curvature (mm)	[83.5in]

# MASK_HT0 is calculated as ZPT_YM * tan (M_ANGLE) - MASK_ZINTCPT (-0.4in)
define	MASK_HT0 3.378		# Height (mm) above datum at SMCS y=0

define	PPLDIST	20018.4D0	# Distance to Telescope Exit Pupil (mm) (REF)

# define	MM_PIX	0.15767		# mm per pixel
# define	ASEC_PIX 0.21		# arcsec per pixel (unbinned)

define	PIX_SZ		0.015		# pixel sz (mm)
define	CCDXPIX		 2048.		# 2048. pixels wide
define	CCDYPIX		 4096.		# 4096. pixels deep
define	CCDXEDG		0.154		# 84 micron to P3 + 70 um to sawline
define	CCDYEDG		0.070		# 40 micron to P3 + 30 um to sawline
define	NOMXGAP		   1.		# 1 mm gap
define	NOMYGAP		 0.1		# 0.1 mm gap
define	FCSYEDG		0.140		# 140 um active to sawline
define	FCSYPIX		  600.		# 600 pixels high as used

define		SZ_INLN	128	# Maximum length of input line
define		SZ_ID	32	# Maximum length of ID string
define		CODE_RF	0	# Priority code for reference (addn) objects
define		CODE_GS	-1	# Priority code for guide star
define		CODE_AS	-2	# Priority code for alignment star

