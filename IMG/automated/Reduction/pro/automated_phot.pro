;+
; NAME:
;   automated_phot
;
; PURPOSE:
; Procedure that performes automated photometry
;
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; 
; 
; OPTIONAL INPUTS:
;  
; 
; OUTPUTS: 
;
; The offset to apply
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
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


PRO automated_phot, night, FLWO=FLWO, filter

  

  ;;Set up paths and files
  cd, getenv('AUTOM_DIR')+'nights/'+night+'/Sci'           
  output = getenv('AUTOM_DIR')+'/lightcurve/all_photometry_'+filter+'.dat' ;;Sets output file
  log = getenv('AUTOM_DIR')+'log/'+night+'_photometry_'+filter+'.log'     ;;Set name of the photometrylog file
  
  openw, lun_out, output, /APPEND, /get_lun 
  openw, lun_log, log, /get_lun           
  
  ;;Store array of fits files found in the Sci directory in variable filesfound
  spawn, 'ls '+'*'+filter+'.fits', filesfound
  errortest = SIZE(filesfound)
  IF (errortest[0] EQ 0) THEN BEGIN
     SPLOG,'Error no images with this filter (or something worse)!'
     RETURN
  ENDIF
  jd=filesfound
  object=filesfound
  
  ;;Make a list of the science files and JD for each file for the current night
  readcol, getenv('AUTOM_DIR')+'/targets/new_old_names.txt', new_names, old_names, $
           COMMENT='#', FORMAT='A,A' ;;Read in new/old blazar names

  for i=0, n_elements(filesfound)-1 do begin
     header = xheadfits(filesfound[i]) ;;Read header from image file
     IF ( N_ELEMENTS(FLWO) EQ 1 ) THEN BEGIN
        jd[i] = FXPAR(header,'MJD')
     ENDIF ELSE BEGIN
        jd[i] = fxpar(header, 'MJD')
     ENDELSE
     object[i] = fxpar(header, 'OBJECT')
     ;;find if old name
     j = WHERE(old_names EQ STRCOMPRESS(object[i], /REMOVE_ALL)) 
     if j ne [-1] then object[i] = new_names[j]
  endfor

  
  ;;Store number of blazars found for this night	
  nobj = n_elements(object)     
  ;;Read in list of template files	
  readcol, getenv('AUTOM_DIR')+'archive/templates.lst', templates, FORMAT='A' 

  
  ;;Loop through each image for the night and perform photometry on it
  ;;Various error catching if statements prevents failure of the pipeline and logs the errors.
  ;;Most errors occur when grab_image_photometry.pro returns a -1
  FOR k=0, nobj-1 DO BEGIN

     ;;Check if template exists, result is -1 if not
     find_template = where(templates eq STRCOMPRESS(object[k],/REMOVE_ALL)+'.fits')
     ;;Read in image to check for wcs and also for making thumbnails
     ;;REMEBER to change the lower case r back to a capital. Make it more general.
     img_hdr = headfits(getenv('AUTOM_DIR')+'nights/'+night+'/Sci/sci_'+STRCOMPRESS(object[k],/REMOVE_ALL)+'_'+filter+'.fits',EXTEN=0) 
     ;;Read header to check if WCS exists
     extast, img_hdr, astr, wcs_test 
     PRINT,wcs_test
     ;;If a WCS fit or template exists then we can perform photometry

     if find_template ne [-1] or wcs_test eq 2 then begin 

        ;;Do automated photometry to get calibrated magnitudes for the current object
        photometry = automated_grabphot(STRCOMPRESS(object[k],/REMOVE_ALL), night, filter) 
        ;;Output results to all_photometry.dat and the night's
        ;;photometry log if successful
       

        if photometry ne [-1] then begin 

           ;;Output results to all_photometry.dat
           printf, lun_out, object[k], jd[k], photometry, $ 
                   format='(A25,F14.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)'

           if photometry[0] gt 70.0 or photometry[1] gt 70.0 or photometry[4] gt 70.0 then begin
              printf, lun_log, object[k]+' photometry failed. No new data plotted.'

           endif else printf, lun_log, object[k]+' photometry successful!'

           ;;Error on bad image

        endif else printf, lun_log, object[k]+' photometry failed.  Unable to find stars in image. No new data plotted.' 

     endif else printf, lun_log, object[k] + ' has no template and no WCS fit.  Unable to perform photometry, create a template! No new data plotted.'
 
  endfor
  
  ;;back to main
  cd, getenv('AUTOM_DIR')
  FREE_LUN,lun_out
  FREE_LUN,lun_log
  
END
