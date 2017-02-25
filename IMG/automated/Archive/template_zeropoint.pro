
;Program calculates zero point for an image
PRO template_zeropoint, obj_name, img_ra, img_dec, img_mag, output_zeropoint, output_zeropoint_stddev, field=field, plot=plot
readcol, '/b/Blazars/papers/DataI/Lists/fields.lst', field_names, sdss, usno, cm14, FORMAT='A,A,A,A', comment='#'

;If else ladder to determine which catalogue we should use
IF NOT KEYWORD_SET(field) THEN BEGIN
	catalogue='usno' ;Default catalogue that  covers all possible obects
	index = WHERE(obj_name EQ field_names) ;Find the index in the fields list for this object
	IF cm14[index] EQ 'Yes' THEN catalogue='cm14' ;2ndary catalogue if there is nothing in the sdss
	IF sdss[index] EQ 'Yes' THEN catalogue='sdss' ;Primary catalogue for getting photometry
ENDIF ELSE catalogue=field
cat_file = '/b/Blazars/papers/DataI/Lists/ref_stars/'+obj_name+'_'+catalogue+'.lst' ;set the refstars file name after choosing the right catalogue automatically
READCOL, cat_file, cat_ra, cat_dec, cat_mag, COMMENT='#', FORMAT='D,D,D' ;Reads in RA, DEC and known magnitudes for the ref stars

n_cat = N_ELEMENTS(cat_ra) ;Number of catalogue stars found
img_stars_in_cat = MAKE_ARRAY(n_cat, /INT) ;Array to store indicies of cat. stars found in image
zp = MAKE_ARRAY(n_cat, /FLOAT) ;Store zero point for each star found in catalogue
zp_stddev = MAKE_ARRAY(n_cat, /FLOAT) ;Store zero point stddev for each star found in catalogue
rad = 0.00010D ;Radius (in degrees) to coorelate the RA and DEC of stars in image to stars in catalogue
FOR i = 0, n_cat-1 DO BEGIN ;Loop coorelates stars from catalogue with stars in the image
	j = WHERE(img_ra GT cat_ra[i] - rad AND img_ra LT cat_ra[i] + rad AND $ ;Search for catalogue star in RA and Dec.
			  img_dec GT cat_dec[i] - rad AND img_dec LT cat_dec[i] + rad, nj)
	IF nj GT 1 THEN BEGIN ;If more than one star is found in the radius, pick the brightest
		holder = MIN(img_mag[j], brightest_star) ;Find index of brightest_star
		j = j[brightest_star] ;Set j to be only the brightest star found in the search radius
		nj = 1 ;Set stars found to 1
	ENDIF
	IF nj EQ 1 THEN  BEGIN
		img_stars_in_cat[i] = 1 ;Show that the catalogue star was found in the image
		zp[i] = img_mag[j] - cat_mag[i] ;Calculate zero piont
		;;zp_stddev[i] = SQRT(img_stddev[j]^2+cat_stddev[i]^2) ;Calculate standard deviation for the zero point
	ENDIF
ENDFOR
cat_mag = cat_mag[WHERE(img_stars_in_cat EQ 1)] ;Only use catalogue stars found in image
cat_stars_in_img = WHERE(zp NE 0) ;Catalogue stars found in image
zp = zp[cat_stars_in_img] ;Remove zeros
PRINT, 'Found '+STRING(N_ELEMENTS(zp))+' catalogue stars in the image!'
PRINT, 'Now time to make some cuts to get the zero piont.'

;;;zp = img_mag[cat_stars_in_img] - cat_mag[img_stars_in_cat]

IF KEYWORD_SET(plot) THEN BEGIN
	;Make plot for making zero point cuts
	plot, cat_mag, zp, psym=1, xtitle='Catalogue R-Mag.', ytitle='Zero Point',$  ;examine relation between the zero point and known magnitudes
			xr=[9.,23.], yr=[-23.2, -19.]
ENDIF

;Try to automatically do mag and zp cuts
DJS_ITERSTAT, cat_mag, SIGREJ=2.5, MASK=mag_mask, MEAN=mag_mean, SIGMA=mag_sigma ;Make a 3-sigma cut
mag_cut = WHERE(mag_mask EQ 1) ;Make cuts
DJS_ITERSTAT, zp[mag_cut], SIGREJ=2.5, MASK=zp_mask, MEAN=zp_mean, SIGMA=zp_sigma
zp_cut = WHERE(zp_mask EQ 1)

zp_mean = AVG(zp[mag_cut[zp_cut]])
zp_stddev = STDDEV(zp[mag_cut[zp_cut]])
IF KEYWORD_SET(plot) THEN BEGIN
	oplot, [MIN(cat_mag[mag_cut[zp_cut]]), MIN(cat_mag[mag_cut[zp_cut]])], [-100,100] ;Plot mag cuts
	oplot, [MAX(cat_mag[mag_cut[zp_cut]]), MAX(cat_mag[mag_cut[zp_cut]])], [-100,100]
	oplot, [-100,100], [MIN(zp[mag_cut[zp_cut]]), MIN(zp[mag_cut[zp_cut]])] ;Plot zp cuts
	oplot, [-100,100], [MAX(zp[mag_cut[zp_cut]]), MAX(zp[mag_cut[zp_cut]])]
	answer = '' ;Initialize
	REPEAT BEGIN
		READ, answer, PROMPT='Do you want to redo the cuts manually? (y/n) '  ;Ask user if they want to redo the cuts manually
		IF answer EQ 'y' THEN BEGIN
			READ, minzp, maxzp, PROMPT='What is the min, max zero points? '
			READ, minmag, maxmag, PROMPT = 'What is the min, max mag.? '
			mag_cut = WHERE(cat_mag GE minmag AND cat_mag LE maxmag)
			zp_cut = WHERE(zp[mag_cut] GE minzp AND zp[mag_cut] LE maxzp)
			oplot, [MIN(cat_mag[mag_cut[zp_cut]]), MIN(cat_mag[mag_cut[zp_cut]])], [-100,100] ;Plot mag cuts
			oplot, [MAX(cat_mag[mag_cut[zp_cut]]), MAX(cat_mag[mag_cut[zp_cut]])], [-100,100]
			oplot, [-100,100], [MIN(zp[mag_cut[zp_cut]]), MIN(zp[mag_cut[zp_cut]])] ;Plot zp cuts
			oplot, [-100,100], [MAX(zp[mag_cut[zp_cut]]), MAX(zp[mag_cut[zp_cut]])]
			zp_mean = AVG(zp[mag_cut[zp_cut]])
			zp_sigma = STDDEV(zp[mag_cut[zp_cut]])
		ENDIF
	ENDREP UNTIL answer EQ 'n'
ENDIF

;SET THE OUTPUT VARIABLES
output_zeropoint = zp_mean
output_zeropoint_stddev = zp_sigma
	
RETURN
END
