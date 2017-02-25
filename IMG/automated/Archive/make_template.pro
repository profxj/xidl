;This program is a test for making a template for a single blazar from a good Super-LOTIS fits image
;from which can later be read from for automating the pipeline (mainly fitting WCS and doing photometry)
PRO make_template, file, FIELD=field

image = XMRDFITS(file, 0) ;Read in image from fits file
header = XHEADFITS(file) ;Read in header from fits file
blazar = STRCOMPRESS(FXPAR(header, 'OBJECT'), /REMOVE_ALL) ;Grab name of blazar from header
jd = FXPAR(header, 'JD')

PRINT, 'Check image to find out what the source aperature and background radii are needed for the blazar.' 
XATV, file, /BLOCK ;Allow user to inspect image
READ, blz_src_aper, blz_inner_bkground, blz_outer_bkground, PROMPT='Enter for blazar (src aper, inner bkground, outer bkground):' 

;Search for all stars that can be found in the image and stores their
;x,y coords, RA & DEC, instrumental magnitude as determined from APER
FIND, image, x, y, flux, sharp, round, 0.5, 10.0, [-1.0,1.0], [0.2,1.0] ;Find all objects in field
gain = FXPAR(header, 'GAIN')
APER, image, x, y, instr_mag, instr_stddev, sky, skyerr, gain ;Perform aperature phot. on all objects in field

;Eliminate the field stars with bad photometry
good_phot = WHERE(instr_mag NE 99.999) ;Find all stars in image where the photometry failed
instr_mag = instr_mag[good_phot]
instr_stddev = instr_stddev[good_phot]
x = x[good_phot]
y = y[good_phot]
sky = sky[good_phot]
skyerr = skyerr[good_phot]
sharp = sharp[good_phot]
round = round[good_phot]

instr_mag = instr_mag - 25 ;Get rid of annoying default zero point of 25 from the instrumental magnitudes

;Check if WCS is fit in image and have the user fit it manually if it is not
extast, header, astr, wcs_test ;Read header to check if WCS exists
IF wcs_test EQ 2 THEN BEGIN
	print, 'WCS found!'
ENDIF ELSE BEGIN
	print, 'WCS not found!  Need to manually set WCS'
	center_ra = FXPAR(header, 'RA') ;Get RA and Dec. 
	center_dec = FXPAR(header, 'DEC')
	radius = 17 ; The FOV radius for Super-LOTIS is 17 arcmin
	REPEAT BEGIN ;Repeat until user is satisfied with WCS fit
		QUICK_ASTRO, file, center_ra, center_dec, radius, SETHEAD=header ;Set header to now include WCS fit
		PRINT, 'Check to make sure the WCS fit is correct'
		XATV, image,HEADER=header, /BLOCK ;Allows user to inspect WCS fit
		answer = ''
		READ, answer, PROMPT='Redo WCS fitting? (y/n) '
	ENDREP UNTIL answer EQ 'n'
ENDELSE
;Now that the WCS is guaranteed to be fit, we will convert the x,y pixel coordinates into ra and dec
XYAD, header, x, y, ra, dec ;Convert x,y for each star into ra,dec


;Get zero point for the image to determine the "photometric" magnitudes for each star
IF KEYWORD_SET(field) THEN template_zeropoint, blazar, ra, dec, instr_mag, zp, zp_stddev, Field=field, /PLOT $
	ELSE template_zeropoint, blazar, ra, dec, instr_mag, zp, zp_stddev, /PLOT
phot_mag = instr_mag - zp ;Calculate photometric magnitudes for each star
phot_stddev = SQRT(zp_stddev^2+instr_stddev^2) ;Calculate the error on the photometric magnitude for each star

;Find reference star in template and store it's magnitude
radius = 0.0014 ;RA and dec to search In degrees (~5 arcsec.)
READCOL, '/b/Blazars/data/target_list.dat', target, blz_ra, blz_dec, str_ra, str_dec, FORMAT='A,A,A,A,A', COMMENT='#'
tar = WHERE(target EQ blazar) ;Find index in target array where it is the same as the name of the blazar we are looking at
X_RADEC, [blz_ra[tar], str_ra[tar]], [blz_dec[tar], str_dec[tar]], phot_ra, phot_dec ;Store blazar and ref. stars' ra and dec in the phot arrays

;Now cross coorelate reference star with template
ref_star_index = WHERE(ABS(ra - phot_ra[1]) LE  radius) ;Find all stars in template within 5 pixels in x direction that match reference star
ref_star_index = ref_star_index[WHERE(ABS(dec[ref_star_index] -phot_dec[1]) LE  radius)] ; same for y direction
IF n_elements(ref_star_index) EQ 1 THEN BEGIN
	PRINT, 'Found the reference star in template!'
ENDIF ELSE BEGIN
	PRINT, 'ERROR: Could not find reference star in template'
ENDELSE

;Store everything as a IDL structure.  It essentially becomes a template Super-LOTIS image from which
;the WCS fitting and photometry can be done automatically!
field = {x:x, y:y, ra:ra, dec:dec, mag:phot_mag, stddev:phot_stddev} ;Create field structure
FXADDPAR, field_hdr, 'Object', blazar, 'Name of object'
FXADDPAR, field_hdr, 'JD', jd, 'Julian date image was taken'
FXADDPAR, field_hdr, 'REF_MAG', phot_mag[ref_star_index[0]], 'Calibrated mag. for ref. star'
FXADDPAR, field_hdr, 'REF_SIG', phot_stddev[ref_star_index[0]], 'Std. dev for ref. star mag'
FXADDPAR, field_hdr, 'ZP', zp, 'Mean zero point of image'
FXADDPAR, field_hdr, 'ZP_SIGMA', zp_stddev, 'Std. dev. for the zero point calcluation'
FXADDPAR, field_hdr, 'FILE', file, 'Original image'
FXADDPAR, field_hdr, 'SRC_APER', blz_src_aper, 'Source aperture for the blazar'
FXADDPAR, field_hdr, 'IN_BK', blz_inner_bkground, 'Inner radius for blz. background'
FXADDPAR, field_hdr, 'OUT_BK', blz_outer_bkground, 'Outer radius for blz. background'

MWRFITS, field, blazar+'.fits', field_hdr, /CREATE ;Store IDL structure
SPAWN, 'ls *.fits > templates.lst' ;Store list of current template files

PRINT, 'zp=', zp
PRINT, 'stddev=', zp_stddev

CLOSE, /ALL
RETURN
END
