;+
; NAME:
;   automated_grabphot
;
; PURPOSE:
; 
; Program performs photometry on a given image.  If image is not WCS
; fit, it tries to fit astrometry using an existing template.
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
; It returns an array with all the photometry results (see end of program for format information)
; If it fails it will return -1 which indicates an error occured (usually a poor image)
;
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
; result = grab_image_photometry('1ES_0806+524', '110222')
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;           2011  Written by KK
;    05-Jun-2012  Revised by MF
;
;-
;------------------------------------------------------------------------------
;; 


function automated_grabphot, object, night, filter, radius=radius, old=old

  
  ;;Reading in data and initialize all the variables
  ;;Uses old name if it is inputted by the variable 'old'
  if keyword_set(old) then name=old[0] else name=object 
  ;;Name of the image to perform photomtery on if filter is R
  ;;changed * to _r
  file = getenv('AUTOM_DIR')+'nights/'+night+'/Sci/sci_'+name+'_'+filter+'.fit*[sz]' 
  
  ;;-----------------------------
  ;;Establish Aperature annulus
  ;;-----------------------------

  image = MRDFITS(file, 0, img_hdr)                                  ;Read in image and header
  IF NOT KEYWORD_SET(radius) THEN radius=5                           ;Set search radius in pixels
  aperature = 15                                                     ;Number of pixels for aperature photometry
  skyannulus = [20,30]                                               ;Pixels describing sky annulus radii
  template = getenv('AUTOM_DIR')+'archive/'+object+'.fits'           ;Set path for template (if it exists)
 

  ;;Load in information from the template into into structure "field" 
  field = MRDFITS(template, 1, fld_hdr) 
  

  ;;If a template (ideally) exists, grab source and background aperatures form the template, 
  ;;else go to default if no template exists
  if keyword_set(fld_hdr) then begin  
     src_aper = FLOAT(FXPAR(fld_hdr, 'SRC_APER'))            ;Grab source apture region radius
     src_inner_bkground = FLOAT(FXPAR(fld_hdr, 'IN_BK'))     ;Grab background region inner radius annulus
     src_outer_bkground = FLOAT(FXPAR(fld_hdr, 'OUT_BK'))    ;Grab background region outer radius annulus
  endif else begin 
     src_aper = 15              ;Default source aperture region radius to 15 pixels
     src_inner_bkground = 20    ;Default  background annulus inner radius to 20 pixels
     src_outer_bkground = 30    ;Default background annulus outer radius to 30 pixels
  endelse
  
  ;;-----------------------------
  ;;Determine good photometry
  ;;-----------------------------
  
  n_brtstr = 40                 ;;Set the number of brightest stars to match
  gain  = 1.                    ;;FLOAT(FXPAR(img_hdr, 'GAIN')) ;already in e'

  ;;Find all objects in the image and perform aperture photometry 
  ;;make a new fits image and put a white dot where each of these x
  ;;and y positions are.
  find, image, x, y, flux, sharp, round, 0.5, 10.0, [-1.0,1.0], [0.2,1.0], /SILENT
  
  aper, image, x, y, img_mag, img_stddev, sky, skyerr, gain, aperature, skyannulus, [-5000,5000], /SILENT
  img_mag = img_mag - 25  ;;Get rid of the annoying default zero point
  ;;Use only stars in image with good photometry
  good_phot = WHERE(img_mag LE 73.999) ;;Finds indicies of only the good photometry
  x = x[good_phot]
  y = y[good_phot]
  img_mag = img_mag[good_phot]
  img_stddev = img_stddev[good_phot]
  imageTest = image
  imageTest2 = image

  
  ;;-------------------------------------------
  ;;Find brightest stars in Template and Image
  ;;--------------------------------------------

;If no WCS fit exists the program will scale and pan (find x,y offsets) the template stars to match the image
;and then use starast.pro to try to fit the astrometry
;The code that follows is long and complicated and only needs running when the astrometry is not 
;initially properly fit, but I will try my best to explain it.
  extast, img_hdr, astr, wcs_test ;Read header to check if WCS exists
  IF wcs_test NE 2 THEN BEGIN 
     ;;Check to make sure n_brtstr bright stars can be found in image & template	
     IF N_ELEMENTS(img_mag) LT n_brtstr OR N_ELEMENTS(field.mag) LT n_brtstr THEN BEGIN

        SPLOG,'not enough stars!'
        RETURN, -1 
     ENDIF
     fld_brtstr = SORT(field.mag)             ;Sort template stars by magnitude
     fld_brtstr = fld_brtstr[0:n_brtstr-1]    ;Find brightest stars in template
     img_brtstr = SORT(img_mag)               ;Sort image stars by magniutde
     img_brtstr = img_brtstr[0:n_brtstr-1]    ;Find brightest stars in image
     
     ;;--------------------------------------------------
     ;;FIND SCALE DIFFERENCES BETWEEN TEMPLATE AND IMAGE
     ;;--------------------------------------------------

     binsize = 5                                               ;Set bin size for the up and coming histograms
     scale_range = [0.90,1.10]                                 ;Range of rotations (in radians)
     scale_step = 0.001                                        ;Step between each test of rotation (in radians)
     n_steps = (scale_range[1] - scale_range[0])/scale_step ;Calculate number of steps for finding the rotation angle
     store_offset_binsizes_x = [0.,0.] ;Initialize 2 element array  that stores the max bin size found in the offsets for each step in the rotation
     store_offset_binsizes_y = [0.,0.]
     FOR k=0, n_steps-1 DO BEGIN                      ;Scale bit by bit and search for best fit
        scale = scale_range[0] + k*scale_step         ;Angle for this step to test
        result_x = field.x[fld_brtstr]*scale
        result_y = field.y[fld_brtstr]*scale
        delta_x = FLTARR(n_brtstr, n_brtstr, /NOZERO) ;Make arrays that stores the the differences in x,y between image stars and template stars
        delta_y = FLTARR(n_brtstr, n_brtstr, /NOZERO)
        FOR i=0, n_brtstr-1 DO BEGIN ;Find delta x,y for the brightest image stars for all the brightest template stars
           delta_x[i,*] = x[img_brtstr[i]] - result_x
           delta_y[i,*] = y[img_brtstr[i]] - result_y
        ENDFOR
        delta_xbins = HISTOGRAM(delta_x, BINSIZE=binsize, LOCATION=delta_xbin_values) ;Bin up delta x y
        delta_ybins = HISTOGRAM(delta_y, BINSIZE=binsize, LOCATION=delta_ybin_values)
                                ;sets the minimum offsets
        IF store_offset_binsizes_x[0] LT MAX(delta_xbins) THEN BEGIN
           store_offset_binsizes_x=[MAX(delta_xbins), MAX(delta_ybins)]
           scale_x = scale
        ENDIF
        IF store_offset_binsizes_y[1] LT MAX(delta_ybins) THEN BEGIN
           store_offset_binsizes_y=[MAX(delta_xbins), MAX(delta_ybins)]
           scale_y = scale
        ENDIF
     ENDFOR

     ;;----------------------------------------------
     ;;Apply the scale angle and find x,y, offsets
     ;;----------------------------------------------

     ;;Now the scale angle is found, apply it to the template and then find the x,y-offsets
     delta_x = FLTARR(n_brtstr, n_brtstr, /NOZERO) ;Make arrays that stores the the differences in x,y between image stars and template stars
     delta_y = FLTARR(n_brtstr, n_brtstr, /NOZERO)
     result_x = field.x[fld_brtstr]*scale_x
     result_y = field.y[fld_brtstr]*scale_y
     FOR i=0, n_brtstr-1 DO BEGIN ;Find delta x,y for the brightest image stars for all the brightest template stars
        delta_x[i,*] = x[img_brtstr[i]] - result_x
        delta_y[i,*] = y[img_brtstr[i]] - result_y
     ENDFOR
     
     ;;Bin up delta x,y because the largest bin contains the true x,y offsets of the image from the template
     delta_xbins = HISTOGRAM(delta_x, BINSIZE=binsize, LOCATION=delta_xbin_values) ;Bin up delta x y
     delta_ybins = HISTOGRAM(delta_y, BINSIZE=binsize, LOCATION=delta_ybin_values)
     max_xbin = WHERE(delta_xbins EQ MAX(delta_xbins)) ;Find the location of the largest bins
     max_ybin = WHERE(delta_ybins EQ MAX(delta_ybins))
     max_xbin_value = delta_xbin_values[max_xbin] ;Get values of largest bins
     max_ybin_value = delta_ybin_values[max_ybin]
     IF N_ELEMENTS(max_xbin) NE 1 THEN max_xbin_value = AVG(max_xbin_value)
     IF N_ELEMENTS(max_ybin) NE 1 THEN max_ybin_value = AVG(max_ybin_value)
     xbin_points = WHERE(delta_x GE max_xbin_value[0] AND delta_x LE max_xbin_value[0]+5.0) ;Store indicies in delta x,y found only in the largets bins
     ybin_points = WHERE(delta_y GE max_ybin_value[0] AND delta_y LE max_ybin_value[0]+5.0)
     x_offset = MEDIAN(delta_x[xbin_points]) ;Take median of largest bins to be the x,y offsets, should be VERY robust!
     y_offset = MEDIAN(delta_y[ybin_points])
     
     ;;Correct for scale and x,y offsets to make matching the stars much easier in the later code
     field.x = field.x*scale_x + x_offset ;Store the new coords. for the template after correcting for scaling and x,y offsets
     field.y = field.y*scale_y + y_offset
     
     ;;Match stars found in image with stars found in template
     nfld = N_ELEMENTS(field.x) ;count number of stars in template
     nimg = N_ELEMENTS(x)       ;count number of stars found in image
     index_in_img = MAKE_ARRAY(nfld, /INTEGER)
     FOR i=0, nfld-1 DO BEGIN   
        ;;Try to coorelate each star from the template with a star
        ;;found in the image
        ;;Find all stars in image withen +/- radius of i'th template star in x-direction
        search_x = WHERE(x GE field.x[i] - radius AND x LE field.x[i] + radius) 
        ;;Find all stars in image withen +/- radius of i'th template star in y-direction
        search_y = WHERE(y GE field.y[i] - radius AND y LE field.y[i] + radius) 
        match = [-1]                             ;Erase and initialize match
        ;;Loop through each piece of search_x and see if it matches anything in search_y
        FOR j=0, N_ELEMENTS(search_x)-1 DO BEGIN 
           find_index = WHERE(search_y EQ search_x[j]) ;Find the index in search_y that is equal to search_x
           IF find_index NE [-1] THEN BEGIN            ;Run only if an index is found
              IF match EQ [-1] THEN match=[search_y[find_index]] $ ;First match found
              ELSE match = [match, search_y[find_index]]           ;Any subsequent matches found
           ENDIF
        ENDFOR
        ;;If one match is found store it in array index_in_img
        IF match NE [-1] AND N_ELEMENTS(match) EQ 1 THEN index_in_img[i] = match 
     ENDFOR
    
     ;;Use only found template stars by removing template stars not found in image
     nozero = WHERE(index_in_img NE 0)
     IF nozero EQ [-1] THEN BEGIN 
        SPLOG,'No stars correspond to the template!'
        RETURN, -1 ;If no stars are found in image that corrispond to the template
     ENDIF
     index_in_img = index_in_img[nozero]
     x = x[index_in_img]        ;Ignore any stars not found in image
     y = y[index_in_img]
     img_mag = img_mag[index_in_img]
     img_stddev = img_stddev[index_in_img]
     field = {x:field.x[nozero], y:field.y[nozero], ra:field.ra[nozero],$
              dec:field.dec[nozero], mag:field.mag[nozero],stddev:field.stddev[nozero]}
     
     ;;Find brighest stars in first three quadrants so the stars are
     ;;likely well seperated when fitting the wcs with starast
     quad1 = WHERE(x LT 1024.0 AND y LT 1024.0)
     quad2 = WHERE(x LT 1024.0 AND y GE 1024.0)
     quad3 = WHERE(x GE 1024.0 AND y LT 1024.0)
     quad4 = WHERE(x GE 1024.0 AND y GE 1024.0)
     IF quad1 EQ [-1] OR quad2 EQ [-1] OR quad3 EQ [-1] OR quad4 EQ [-1] THEN BEGIN

        SPLOG,'Error catching!'
        RETURN, -1              ;Error catching

     ENDIF

     radii = SQRT((x-1024)^2 + (y-1024)^2) ;Calculate radius from image center for all stars
     ;;Find the star with the max radius from the center in a given quadrant
     hldr = MAX(radii[quad1], maxrad1)     
     hldr = MAX(radii[quad2], maxrad2)
     hldr = MAX(radii[quad3], maxrad3)
     hldr = MAX(radii[quad4], maxrad4)
     wcs_fitting_stars = [quad1[maxrad1[0]], quad2[maxrad2[0]], quad3[maxrad3[0]]]
     
     ;;Fit WCS (if it dozesn't already exist)
     STARAST, field.ra[wcs_fitting_stars], field.dec[wcs_fitting_stars],$
              x[wcs_fitting_stars], y[wcs_fitting_stars], HDR = img_hdr
  ENDIF

  

   
  ;;Grab zero point from image using the stars in the image and comparing them to stars from a catalogue
  xyad, img_hdr, x, y, img_ra, img_dec ;;Grab RA and Dec. for all the stars found in the image
  automated_imagezp, object, img_ra, img_dec, img_mag, img_stddev, img_zp, img_zp_stddev, zp_numused 
  if not keyword_set(img_zp) then  BEGIN 
     SPLOG,'cant find zero point'
     
     return, -1 ;;Catch error when unable to find zero point
  ENDIF

  ;;Read in coordinates for blazar and reference star in order to do photometry
  READCOL, getenv('AUTOM_DIR')+'/targets/target_list.dat', target, blz_ra, blz_dec, str_ra, str_dec, $
           FORMAT='A,A,A,A,A', COMMENT='#'
  ;;Find index in target array where it is the same as the name of the blazar we are looking at
  tar = WHERE(target EQ object) 
  ;;Store blazar and ref. stars' ra and dec in the phot arrays
  X_RADEC, [blz_ra[tar], str_ra[tar]], [blz_dec[tar], str_dec[tar]], phot_ra, phot_dec 
  
  ;;Now grab predicted pixel locations for blazar and reference star and perform relative photometry on them
  source_rad = 15                                                ;Set radius for blazar/ refstar 
  ;;Convert RA and DEC for blz and ref star into predicted x,y coods. in img
  ADXY, img_hdr, phot_ra, phot_dec, phot_approx_x, phot_approx_y 
  ;;Find the centroid of the blazar and ref. star from the predicted x,y coords.
  CNTRD, image, phot_approx_x, phot_approx_y, phot_cntrd_x, phot_cntrd_y, source_rad, /SILENT 
  APER, image, phot_cntrd_x[0], phot_cntrd_y[0], blazar_mag, blazar_stddev, sky, skyerr, gain, src_aper, $
        [src_inner_bkground, src_outer_bkground], [-5000,5000], /SILENT ;Photometry of ref. star
  APER, image, phot_cntrd_x[1], phot_cntrd_y[1], refstar_mag, refstar_stddev, sky, skyerr, gain,$
        aperature, skyannulus, [-5000,5000], /SILENT ;Photometry of ref. star
  phot_mag = [blazar_mag, refstar_mag]
  phot_stddev = [blazar_stddev, refstar_stddev]
  phot_mag = phot_mag - 25      ;Get rid of the default 25 zero point added to the photometry by APER
  
 ;;;Now cross coorelate refernece star with template
  ;;ref_pos = WHERE(ABS(field.x - phot_cntrd_x[1]) LE  radius) ;Find all stars in template within 5 pixels in x direction that match reference star
  ;;ref_pos = ref_pos[WHERE(ABS(field.y[ref_pos] - phot_cntrd_y[1]) LE  radius)] ; same for y direction
  ;;IF n_elements(ref_pos) EQ 1 THEN BEGIN
  ;;	PRINT, 'Found the reference star in template!'
  ;;ENDIF ELSE BEGIN
  ;;	PRINT, 'ERROR: Could not find reference star in template
  ;;	return, -1
  ;;ENDELSE

  ;;Now we calibrate the blazar's magnitude
  ;;calib_refstar_mag = FLOAT(FXPAR(fld_hdr, 'REF_MAG'))
  ;;calib_refstar_stddev = FLOAT(FXPAR(fld_hdr, 'REF_SIG'))
  ;;calib_blz_mag = phot_mag[0] - phot_mag[1] + calib_refstar_mag
  ;;calib_blz_stddev = SQRT(phot_stddev[0]^2 + phot_stddev[1]^2 + calib_refstar_stddev^2)
  
  calib_refstar_mag = phot_mag[1] - img_zp
  calib_refstar_stddev = SQRT(phot_stddev[1]^2 + img_zp_stddev^2)
  calib_blz_mag = phot_mag[0] - img_zp
  calib_blz_stddev = SQRT(phot_stddev[0]^2 + img_zp_stddev^2)
  

  ;;If photometry is successful then output a thumbnail of the image
  make_thumbnail, image, getenv('AUTOM_DIR')+'/trigger_webpage/thumbnails/'+object+'.jpg' ;Create thumbnail

  ;;Format is
  ;;Calibrated Blazar mag.
  ;;Calibrated ref star mag.
  ;;Calibrated blazar std-dev
  ;;Calibrated ref star std-dev
  ;;Relative Blazar mag
  ;;Relative ref star mag
  ;;Relative blazar std-dev
  ;;Relative ref star std-dev
  ;;Zero point
  ;;Zero point std-dev
  ;;Number of stars used to determine the zero point


SPLOG,'return!'
  RETURN, [calib_blz_mag, calib_refstar_mag, calib_blz_stddev, calib_refstar_stddev, $
           phot_mag, phot_stddev, img_zp, img_zp_stddev, zp_numused] ;Return resulting photometry


END
