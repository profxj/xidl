define		NPARAM	42		# no of param (last + 24)
define		COL_DST $1[1]		# Collimator tilt error
define		COL_ERR $1[2]		# Collimator tilt error
define		COL_PHI $1[3]		# Phi of Collimator tilt
define		TNT_ANG	$1[4]		# tent mirror angle
define		TNT_PHI	$1[5]		# tent mirror rotation about z
define		GRLINES	$1[6]		# lines per micron
define		ORDER	$1[7]		# order
define		MU	$1[8]		# Grating angle (true)
define		GR_YERR	$1[9]		# Grating roll error (rotation about y)
define		GR_ZERR	$1[10]		# Grating yaw error (rotation about z)
define		CAM_FOC	$1[11]		# camera focus (mm)
define		CAM_ANG $1[12]		# camera angle (2.33 deg)
define		CAM_PHI $1[13]		# camera phi
define		C2	$1[14]		# camera distortion coeff
define		C4	$1[15]		# camera distortion coeff
define		X_OPT	$1[16]		# x position of optical axis (ICS)
define		Y_OPT	$1[17]		# y position of optical axis (ICS)
define		MOS_ROT	$1[18]		# CCD rotation wrt mask (rad)

define		CN_XERR	$1[16+$2*3]	# x-adjustment in pixels
define		CN_YERR	$1[17+$2*3]	# y-adjustment in pixels
define		CN_RERR	$1[18+$2*3]	# rotation error

define	PIX_SZ	0.015			# pixel size in mm
define	NCCD	8

define	ON_CHIP		0
define	OFF_CHIP	1
