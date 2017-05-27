# define	FL_TEL	150393.95D0	# Foc Len of KII WITH 3" PULLBACK (T8_im.zmx)
# THIS FILE NEEDS UPDATING
# NB redfine FL_TEL/TEL_FOC; etc, plate_scales


# define		TEL_FOC	149583.D0		# focal length (mm)
define		R_IMSURF 2133.6D0	# radius of image surface (mm) 
					# NB - should match ray trace program!!
#define		DIST_C0	4.535218e-4 	# distortion term in r (Sutin 3" Ray Tr)
define		DIST_C0	0.0e-4 		# distortion term in r scale should be 0
define		DIST_C2	-1.111311e-8	# distortion term in r**3  (Sutin 3" RT)

define		OBS_LAT	19.8		# Keck latitude
#define		ASECPMM	1.378		# Arcsec/mm in focal plane  [XXX?]
#define		MM_ARCS	0.7252		# mm/arcsec in focal plane NB: TEST!!!
define		MM_ARCS	0.7253		# mm/arcsec in focal plane  (gnomonic)
define		GEN_ZMX	75.0		# General max zenith angle

define		SRV_AZ1	185.0		# Azimuth of Service tower exclusion OLD
define		SRV_AZ2	332.0		# Azimuth of Service tower exclusion OLD
define		SRV_ZMX	53.2		# Z-angle of Service tower exclusion OLD
#define		SRV_AZ1	1.6		# Azimuth of Service tower exclusion OLD
#define		SRV_AZ2	151.0		# Azimuth of Service tower exclusion OLD
#define		SRV_ZMX	54.		# Z-angle of Service tower exclusion OLD

define		LM_OK	0.		# Service platform conflict
define		LM_SRVZ	1.		# Service platform conflict
define		LM_GENZ	2.		# Lower shutter conflict

define		PRESSURE	480.	# Standard Air pressure in mm(Hg) (6xxmb)
define		TEMP		0.	# Standard temperature (C)
