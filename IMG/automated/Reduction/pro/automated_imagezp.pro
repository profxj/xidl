;+
; NAME:
;   automated_imagezp
;
; PURPOSE:
; 
; Compute the zp of an image
;
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; 
; OPTIONAL INPUTS:
;  
; 
; OUTPUTS: 
;
;
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
; 
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;           2011  Written by KK
;       Aug-2012  Revised by MF
;
;-
;------------------------------------------------------------------------------
;; 



PRO automated_imagezp, obj_name, img_ra, img_dec, img_mag, img_stddev, output_zeropoint, $
                       output_zeropoint_stddev, output_numused, field=field, plot=plot

  readcol, getenv('AUTOM_DIR')+'targets/fields.lst', field_names, sdss, usno, cm14,$
           FORMAT='A,A,A,A', comment='#'

  ;;If else ladder to determine which catalogue we should use
  IF NOT KEYWORD_SET(field) THEN BEGIN
     catalogue='usno'                       ;Default catalogue that  covers all possible obects
     index = WHERE(obj_name EQ field_names) ;Find the index in the fields list for this object
     IF index EQ [-1] THEN RETURN ;Ignore images taken of blazars that we have long removed from the pipeline
     IF cm14[index] EQ 'Yes' THEN catalogue='cm14'    ;2ndary catalogue if there is nothing in the sdss
     IF sdss[index] EQ 'Yes' THEN catalogue='sdss'    ;Primary catalogue for getting photometry
  ENDIF ELSE catalogue=field
  ;;set the refstars file name after choosing the right catalogue automatically
  cat_file = getenv('AUTOM_DIR')+'/targets/ref_stars/'+obj_name+'_'+catalogue+'.lst' 
  ;;Reads in RA, DEC and known magnitudes for the ref stars
  READCOL, cat_file, cat_ra, cat_dec, cat_mag, COMMENT='#', FORMAT='D,D,D' 

  ;;Now we want to cross coorelate the stars from the catalogue with stars in the image
  n_cat = N_ELEMENTS(cat_ra)               ;Number of catalogue stars found
  img_stars_in_cat = MAKE_ARRAY(n_cat, /INT) ;Array to store indicies of cat. stars found in image
  zp = MAKE_ARRAY(n_cat, /FLOAT)             ;Store zero point for each star found in catalogue
  zp_stddev = MAKE_ARRAY(n_cat, /FLOAT)      ;Store zero point stddev for each star found in catalogue
  rad = 0.0010D              ;;Radius (in degrees) to coorelate the RA and DEC of stars in image to catalogue
  ;;Loop coorelates stars from catalogue with stars in the image
  FOR i = 0, n_cat-1 DO BEGIN   
     ;;Search for catalogue star in RA and Dec.
     j = WHERE(img_ra GT cat_ra[i] - rad AND img_ra LT cat_ra[i] + rad AND $ 
               img_dec GT cat_dec[i] - rad AND img_dec LT cat_dec[i] + rad, nj)
     ;;If more than one star is found in the radius, pick the brightest
     IF nj GT 1 THEN BEGIN      
        holder = MIN(img_mag[j], brightest_star) ;Find index of brightest_star
        j = j[brightest_star]                    ;Set j to be only the brightest star found in the search radius
        nj = 1                                   ;Set stars found to 1
     ENDIF
     IF nj EQ 1 THEN  BEGIN
        img_stars_in_cat[i] = 1         ;Show that the catalogue star was found in the image
        zp[i] = img_mag[j] - cat_mag[i] ;Calculate zero piont
        ;;zp_stddev[i] = SQRT(img_stddev[j]^2+cat_stddev[i]^2) ;Calculate standard deviation for the zero point
     ENDIF
  ENDFOR
  cat_stars_in_img = WHERE(zp NE 0) ;Catalogue stars found in image
  ;;Error catching for when no stars can be matched between the catalogue and image
  IF cat_stars_in_img EQ [-1] THEN return 
  zp = zp[cat_stars_in_img]     ;Remove zeros
  zp_stddev = zp_stddev[cat_stars_in_img]
  cat_mag = cat_mag[WHERE(img_stars_in_cat EQ 1)] ;Only use catalogue stars found in image
  
  PRINT, 'Found '+STRING(N_ELEMENTS(zp))+' catalogue stars in the image!'
  PRINT, 'Now time to make some cuts to get the zero piont.'
  
  ;;Make plot for making zero point cuts
  IF KEYWORD_SET(plot) THEN BEGIN
     ;;examine relation between the zero point and known magnitudes
     plot, cat_mag, zp, psym=1, xtitle='Catalogue R-Mag.', ytitle='Zero Point',$ 
           xr=[9.,23.], yr=[-23.2, -19.]
  ENDIF
  
  ;;Do some default cuts to get rid of annoying outliers
  ;;Usually the zero point it between -20 and -23
  i = WHERE(zp GE -30.0 AND zp LE -10.0)
  IF i EQ [-1] THEN RETURN      ;Catch error in case no catalouge stars are found
  cat_mag = cat_mag[i]
  zp = zp[i]
  j = WHERE(cat_mag LE 16.0)    ;Use only bright stars to minimize error / systematic effects
  IF j EQ [-1] THEN RETURN      ;Catch error in case no catalouge stars are found
  cat_mag = cat_mag[j]
  zp = zp[j]
  
  ;;Try to automatically do mag and zp cuts
  DJS_ITERSTAT, cat_mag, SIGREJ=2.0, MASK=mag_mask, MEAN=mag_mean, SIGMA=mag_sigma ;Make a 2-sigma cut
  mag_cut = WHERE(mag_mask EQ 1)                                                   ;Make cuts
  IF mag_cut EQ [-1] THEN RETURN
  DJS_ITERSTAT, zp[mag_cut], SIGREJ=2.0, MASK=zp_mask, MEAN=zp_mean, SIGMA=zp_sigma
  zp_cut = WHERE(zp_mask EQ 1)
  
  ;;Calculated  mean and stddev of zero point
  zp_mean = MEAN(zp[mag_cut[zp_cut]])
  zp_sigma = STDDEV(zp[mag_cut[zp_cut]])
  IF KEYWORD_SET(plot) THEN BEGIN
     oplot, [MIN(cat_mag[mag_cut[zp_cut]]), MIN(cat_mag[mag_cut[zp_cut]])], [-100,100] ;Plot mag cuts
     oplot, [MAX(cat_mag[mag_cut[zp_cut]]), MAX(cat_mag[mag_cut[zp_cut]])], [-100,100]
     oplot, [-100,100], [MIN(zp[mag_cut[zp_cut]]), MIN(zp[mag_cut[zp_cut]])] ;Plot zp cuts
     oplot, [-100,100], [MAX(zp[mag_cut[zp_cut]]), MAX(zp[mag_cut[zp_cut]])]
  ENDIF
  
  ;;SET THE OUTPUT VARIABLES
  output_zeropoint = zp_mean
  output_zeropoint_stddev = zp_sigma
  output_numused = N_ELEMENTS(zp[mag_cut[zp_cut]])
  
  RETURN
END
