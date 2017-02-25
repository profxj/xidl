;+
; NAME:
;   automated_replot
;
; PURPOSE:
; Procedure that updates the Blazar Monitoring Program's webpage
;
; CALLING SEQUENCE:
;
;   automated_replot, /UPDATE_WEBPAGE, /EMAIL_ALERT, $
;   /FLUX
;
; INPUTS:
;   
; 
; OPTIONAL INPUTS:
; 
;   /UPDATE_WEBPAGE -- secure copies plots and data files to remote
;                      website
;   /EMAIL_ALERT    -- sends out automated alert email to individual(s)
;   /FLUX           -- Converts SDSS r-mags into flux density
;
; OUTPUTS: 
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
;   READCOL, WRITECOL, PLOT, PLOTERROR, OPLOT, OPLOTERROR, 
;   M_PSOPEN, M_PSCLOSE, MAKE_LIGHTCURVE_PAGE
;
; REVISION HISTORY:
;    05-Jun-2012   Revised by MF
;
;-
;------------------------------------------------------------------------------
;; 


;;Function converts SDSS r-mags into flux density [Jy]
FUNCTION mag_to_flux, m
  b = 1.2e-10                                          ; Softening coefficient for r band
  s = 3631 * 2 * b * SINH( (-ALOG(10)*m/2.5)-ALOG(b) )       ;Convert sinh mag used in SDSS to Jy
  return, s
END

;;Function converts  flux density [Jy] into SDSS r-mags
FUNCTION flux_to_mag, s
  b = 1.2e-10                   ; Softening coefficient for r band
  f_over_f0 = s / 3631          ;Convert Jy to f/f0
  m = -(2.5/ALOG(10)) * ( ASINH(f_over_f0 / (2*b)) + ALOG(b) ) ;Convert f/f0 into magnitude
  return, m
END

;;prints median magnitude according the light curve data
FUNCTION LIGHT_CURVE_MEDIAN_MAG, blz_name
  exists = FILE_TEST(getenv('AUTOM_DIR')+'/data/results/'+blz_name+'_calibrated.dat') 
  IF (exists EQ 1) THEN BEGIN
     READCOL, getenv('AUTOM_DIR')+'/data/results/'+blz_name+'_calibrated.dat', F='X,F,X',mag
     med_mag = MEDIAN(mag)
     IF NOT KEYWORD_SET(med_mag) THEN BEGIN 
        RETURN, -1 
     ENDIF ELSE BEGIN
        RETURN, med_mag 
     ENDELSE
  ENDIF ELSE BEGIN
     RETURN, -1
  ENDELSE
END

;;returns the index of the proper blazar from READCOL arrays to print 
;;their subclass and log nu_syn   
;;sort of redundant function since median is already calculated
;;FUTURE ISSUES
;;Doesn't take into account older names in queue (which is the point of return, -1
;;Edit: add check with old names
FUNCTION subcl_and_lognu, blz_name,names_array
  FOR i=0, N_ELEMENTS(names_array)-1 DO BEGIN
     IF (STRCMP(blz_name,names_array[i],/FOLD_Case) EQ 1) THEN RETURN, i
  ENDFOR
  ;;in case no match at all
  RETURN, -1
END


PRO automated_replot, FLUX=flux, UPDATE_WEBPAGE=update_webpage, $
                      EMAIL_ALERT=email_alert,check_outliers=check_outliers
  
  ;;set to false (0) initially if no email wants to be set, else set to true (1)
  alert=0                
  
  
  ;;-----------------------
  ;; Read in data
  ;;-----------------------
  
  ;;read target list
  READCOL, getenv('AUTOM_DIR')+'/targets/target_list.dat', target, ra,  dec,  $
           starra,  stardec, FORMAT='A,A,A,A,A',SKIPLINE=1,/SILENT
  ntarget=N_ELEMENTS(target)
  READCOL, getenv('AUTOM_DIR')+'/targets/new_old_names.txt', new_names, old_names, $
           COMMENT='#', FORMAT='A,A' ;Read in new/old blazar names
  READCOL, getenv('AUTOM_DIR')+'/lightcurve/all_photometry.dat', r_targ, r_jd, r_blaz_abs, r_star_abs,$
           r_blaz_errabs,r_star_errabs,r_blaz_rel,r_star_rel,$
           r_blaz_errrel, r_star_errrel, r_zp ,r_zp_err ,r_zp_cnt,$
           FORMAT='A,D,F,F,F,F,F,F,F,F,F,F,I', /SILENT, COMMENT='#'
  ;;Take negative of zero point so that it correctly matches the standard way to do it in astronomy
  r_zp = -r_zp   
  READCOL, getenv('AUTOM_DIR')+'/trigger_webpage/trigger_database.dat', td_obj, td_date, $
           td_mjd, td_state, td_mag, $
           td_mag_stddev, FORMAT='A,A,D,A,F,F' ;Read in database of triggers 
  READCOL, getenv('AUTOM_DIR')+'/trigger_webpage/triggers.dat', trig_obj, trig_type, $
           trig_unit, trig_dim, trig_bright, $
           FORMAT='A,A,A,F,F', COMMENT='#' ;Read in data file on triggers
  READCOL, getenv('AUTOM_DIR')+'/trigger_webpage/spectra_mjd.dat', specmjd_obj, specmjd_mjd, FORMAT='A,F', $
           COMMENT='#'          ;Read in data file of spectra and their mjd
  READCOL, getenv('AUTOM_DIR')+'/trigger_webpage/subcl_and_lognu.txt', $
           FORMAT='A,A,A',names_array,subclass,lognu
  
  ;;-------------------------------
  ;; Create data & webpage files 
  ;;-------------------------------
  
  ;;Open file that will store most recent trigger results
  OPENW, 3, getenv('AUTOM_DIR')+'/trigger_webpage/recent_triggers.txt' 
  ;;Open webpage that will show all triggers
  OPENW, 4, getenv('AUTOM_DIR')+'/trigger_webpage/recent_triggers.html' 
  ;;Open webpage that will automatically update for test-lightcurves
  OPENW, 5, getenv('AUTOM_DIR')+'/trigger_webpage/All-Light-Curves.html' 
  OPENW, 12, getenv('AUTOM_DIR')+'/trigger_webpage/Light-Curves.html'
  ;;Stores the current median magnitude of each blazar
  OPENW, 10, getenv('AUTOM_DIR')+'/trigger_webpage/median_mag.txt' 
  ;;OPENW, 11, '/b/Blazars/data/recent_observations.html'
  SPAWN, "date -d -today +'%y%m%d'", today
  
  ;;-------------------------------
  ;; Write to data & webpage files
  ;;-------------------------------
  
  IF KEYWORD_SET(check_outliers) THEN BEGIN
     OPENW, 6, getenv('AUTOM_DIR')+'trigger_webpage/all_outliers_'+STRTRIM(today[0],2)+'.txt'
     printf, 6, '#MJD    DATE    Magnitude   Error   Significance    Moon_RA     Moon_Dec    Phase   Star Rel. Mag   ZP'
  ENDIF
  PRINTF, 10, FORMAT='("#Blazar", 10X, "Median", 5X, "High Thres.", 5X, "Low Thres.")'
  PRINTF, 5, '<table border = "1">'
  
  ;;CD, '/b/Blazars/data'
  ;;Convert all JD to MJD
  r_jd = r_jd - 2400000.5
  
  SPLOG, 'Done reading in differential photometry, begin calibrations'
  
  ;;Update top of blazar light curve page
  time_since_last_image = CURRENT_MJD()-MAX(r_jd) ;Calculate number of days since last image
  PRINTF, 5, 'It has been '+STRING(time_since_last_image, FORMAT='(F5.1)')+$
          ' days since a blazar was imaged.<br><br>'
  
  ;; Light Curve main page "Light-Curves.html"
  PRINTF, 12, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">'
  PRINTF, 12, '<html xmlns="http://www.w3.org/1999/xhtml">'
  PRINTF, 12, '<head>'
  PRINTF, 12, '<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />'
  PRINTF, 12, '<title>Light Curves</title>'
  PRINTF, 12, '<link href="main.css" rel="stylesheet" type="text/css" />'
  PRINTF, 12, '<meta name="Description" content="" />'
  PRINTF, 12, '</head>'
  PRINTF, 12, '<body>'
  PRINTF, 12, '<div class="header">'
  PRINTF, 12, '<div class="logo">'
  PRINTF, 12, '<h1>Light Curves</h1>'
  PRINTF, 12, '</div>'
  PRINTF, 12, '<div class="navigation">'
  PRINTF, 12, '<ul>'
  PRINTF, 12, '<li><a href="All-Light-Curves.html">All Light Curves</a> </li>'
  PRINTF, 12,'<li><a href="http://scipp.ucsc.edu/~afurniss/private/BlazarMonitoring/thumbnail_page.html">Thumbnails from Most Recent Night</a> </li>'
  PRINTF, 12, '</ul>'
  PRINTF, 12, '</div>'
  PRINTF, 12, '</div>'
  PRINTF, 12, 'It has been '+STRING(time_since_last_image, FORMAT='(F5.1)')+' days since a blazar was imaged.<br><br>'
  
  
  ;; Triggers webpage
  PRINTF, 4, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">'
  PRINTF, 4, '<html xmlns="http://www.w3.org/1999/xhtml">'
  PRINTF, 4, '<head>'
  PRINTF, 4, '<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />'
  PRINTF, 4, '<title>Recent Triggers</title>'
  PRINTF, 4, '<link href="http://scipp.ucsc.edu/~afurniss/private/BlazarMonitoring/main.css" rel="stylesheet" type="text/css" />'
  PRINTF, 4, '<meta name="Description" content="" />'
  PRINTF, 4, '</head>'
  PRINTF, 4, '<body>'
  PRINTF, 4, '<div class="header">'
  PRINTF, 4, '<div class="logo">'
  PRINTF, 4, '<h1>Recent Triggers</h1>'
  PRINTF, 4, '</div>'
  PRINTF, 4, '<div class="navigation">'
  PRINTF, 4, '<ul>'
  PRINTF, 4, '<li><a href="http://scipp.ucsc.edu/~afurniss/private/BlazarMonitoring/Light-Curves.html">Light Curves</a> </li>'
  PRINTF, 4, '</ul>'
  PRINTF, 4, '</div>'
  PRINTF, 4, '</div>'
  PRINTF, 4, 'It has been '+STRING(time_since_last_image, FORMAT='(F5.1)')+' days since a blazar was imaged.<br><br>'
  
  
  
  ;;-------------------------------
  ;; Email Alert
  ;;-------------------------------
  
  ;;If time since last image is over 2 days then email everyone and annoy them about it
  IF time_since_last_image GT 2.5 AND KEYWORD_SET(email_alert) THEN BEGIN
     SPAWN, 'echo "It has been'+STRING(time_since_last_image, FORMAT='(F5.1)')+' days since blazars were imaged by Super-LOTIS. '+$
            'This is an automated alert sent by the pipeline if no blazars have been imaged in over 2 days." | mail -s '+$
            '"Blazars not imaged by Super-LOTIS for '+ STRING(time_since_last_image, FORMAT='(F5.1)')+' days" '+$
            'blazars@ucolick.org'
  ENDIF
  
  
  
  ;;-------------------------------
  ;; Generate plots & data files
  ;;-------------------------------
  
  FOR  tar=0, ntarget-1 DO BEGIN ;Loop through each target blazar
     
     ;;Set the old names to new names in the relative photometry
     i=WHERE(new_names EQ target[tar]) ;Find index for the names list for the current target
     
     ;;If the file all_photometry.dat contains two names for the same object, 
     ;;then pick the "new" one according to new_old_names.txt
     new_name=''
     old_name=''
     IF i NE -1 THEN BEGIN
        ;;Because the photometry file stores incorrect names, having this variable stores the bad name
        new_name=new_names[i[0]] 
        old_name=old_names[i[0]]
        ;;Finds indicies in r_targ wheremthe old name is used
        j=WHERE(r_targ EQ old_names[i[0]]) 
        ;;Sets old names found in r_targ to the new names so they are not confused with each other
        IF j NE [-1] THEN r_targ[j] = target[tar] 
    ENDIF	
     
    
     ;;-------------------------------
     ;;Determine good photometry
     ;;-------------------------------
    
     ;;Find targets from the relative photometry
     r_mytarg = WHERE(STRTRIM(r_targ,2) EQ STRTRIM(target[tar],2),r_numfound) 
     
     ;;Run all code only more than one observation is found
     IF r_mytarg NE [-1] AND N_ELEMENTS(r_mytarg) GE 2 THEN BEGIN 
        ;;Do cuts to the photometry
        ;;Calculate difference in zp vs. median  for cutting out
        ;;nights with the I filter
        ;;Eliminates bad photometry
        delta_zp = r_zp[r_mytarg] - MEDIAN(r_zp[r_mytarg]) 
        goodphot = WHERE(r_blaz_abs[r_mytarg] NE 0. AND ABS(r_blaz_abs[r_mytarg]) LE 40. $
                         AND r_blaz_errabs[r_mytarg] LE 0.45 $
                         AND r_jd[r_mytarg] GT 0. AND delta_zp LT 0.75 AND delta_zp GT -1.0) 
        
        ;;Eliminates relative bad photometry
        relative_goodphot = WHERE(r_star_abs[r_mytarg] NE 0. AND r_star_abs[r_mytarg] LE 40. $
                                  AND r_star_errabs[r_mytarg] LE 0.45 $
        AND r_blaz_abs[r_mytarg] NE 0. AND ABS(r_blaz_abs[r_mytarg]) LE 40. AND $
                                  r_blaz_errabs[r_mytarg] LE 0.45 $
                                  AND r_jd[r_mytarg] GT 0. AND delta_zp LT 0.75 AND delta_zp GT -1.0) 
        
        IF goodphot NE [-1] THEN BEGIN
           r_relative_mytarg = r_mytarg[relative_goodphot]
           r_mytarg = r_mytarg[goodphot] ;eliminate zeros
           
           ;;Find index of current blazar in the trigger arrays	
           t = WHERE(trig_obj EQ target[tar]) 
           ;;Convert fixed triggers in flux to magnitudes
           
           IF t NE [-1] THEN  BEGIN
              IF trig_unit[t] EQ 'flux' AND trig_type[t] EQ 'fixed' THEN BEGIN 
                 trig_dim[t] = flux_to_mag(trig_dim[t])
                 trig_bright[t] = flux_to_mag(trig_bright[t])
              ENDIF
           ENDIF
           
           
           ;;-------------------------------
           ;; Determine trigger level
           ;;-------------------------------

           ;;med_mag is later assigned as the median mag for other plots
           med_mag=median(r_blaz_abs[r_mytarg]) ;Calculate the median magnitude of all the points
           IF t NE [-1] THEN BEGIN
              ;;If trigger type is relative then use +/- from median for tirggers
              IF trig_type[t] EQ 'relative' THEN BEGIN 
                 p1mag = med_mag + trig_dim[t]
                 m1mag = med_mag - trig_bright[t]
              ENDIF ELSE BEGIN  
                 ;;If trigger type is fixed then used the fixed trigger levels, ignoring the median
                 p1mag = trig_dim[t] 
                 m1mag = trig_bright[t]
                 calibrated_median_mag = med_mag ;For subtracting on the relative light curves
              ENDELSE
           ENDIF ELSE BEGIN
              p1mag = med_mag + 1.
              m1mag = med_mag - 1.
              calibrated_median_mag = med_mag
           ENDELSE
           
                                ;Print out any triggers 
                                ;PRINTF, 4, '<b>Most recent blazar triggers from Super-LOTIS</b><br><br>'
           
           end_of_list = n_elements(r_blaz_abs[r_mytarg])-1 ;index of the last measurement taken
           current_mag =  r_blaz_abs[r_mytarg[end_of_list]] ;Finds mag. of most current measurement
           current_sigma = r_blaz_errabs(r_mytarg[end_of_list])
           mag_from_median = current_mag - med_mag
           
           current_jd = max(r_jd[r_mytarg])
           CALDAT, MAX(r_jd) + 2400000.5 +1, ut_month, ut_day, ut_year ;Convert JD to UT date for the next day
           ut_year = STRCOMPRESS(ut_year, /REMOVE_ALL)                 ;Convert to strings with no extra spaces
           ut_month = STRCOMPRESS(ut_month, /REMOVE_ALL)
           ut_day = STRCOMPRESS(ut_day, /REMOVE_ALL)

    ;--------------------------------
    ; Produce possible trigger alert
    ;--------------------------------

           IF current_mag GE p1mag OR current_mag LE m1mag then begin
              
              printf, 4, '<p>'
              printf, 4,  'Possible trigger on: '+repchr(strtrim(target[tar],2),'_',' ')
              printf, 4, 'The MJD is: ' + strcompress(current_jd, /remove_all)+'. '
              days_ago = CURRENT_MJD() - current_jd
              printf, 4, 'The trigger occured '+STRING(days_ago, FORMAT='(F7.1)')+' days ago.'
              printf, 4, '<a href=http://scipp.ucsc.edu/~afurniss/private/BlazarMonitoring/lightcurve_webpage/'+$
                      strtrim(target[tar],2)+'.html>Preview light curve</a></p>'
              
              IF current_mag GE p1mag THEN BEGIN
                 ;;how many sigmas rise above the high trigger threshold
                 significance = abs((p1mag - current_mag) / current_sigma) 
              ENDIF ELSE BEGIN
                 ;;how many sigmas fell below the low trigger threshold
                 significance = abs((current_mag - m1mag) / current_sigma) 
              ENDELSE
              
              IF keyword_set(email_alert) then alert=1 ;set the alert on
           endif
           

           ;;Look for index of any previous confirmed triggers (for
           ;;outputting columns on the light-curves) from the trigger database
           prev_trig = WHERE(td_obj EQ strtrim(target[tar],2), n_prev_trig)
           

           ;;Now make calibrated light curves
           SPLOG, 'Making calibrated light-curve for ',+target[tar]
           x_psopen, getenv('AUTOM_DIR')+'/plots/'+strtrim(target[tar],2)+'_variab_calib.ps', /maxs
           a = WHERE(r_star_abs[r_mytarg] LE 40.0) ;Select only good ref. star photometry
           
           
           ;;--------------------------------------
           ;; dim and bright upper bounds of plots
           ;;--------------------------------------
           
           ;;dim=0.0,
           IF t NE -1 AND trig_dim[t] LT 40.0 THEN BEGIN
              IF STRCMP(STRTRIM(trig_type[t],2),'fixed') EQ 1 AND STRCMP(STRTRIM(trig_unit[t],2),'mag') EQ 1 THEN BEGIN
                 dim=trig_dim[t]
              ENDIF ELSE IF STRCMP(STRTRIM(trig_unit[t],2),'mag') EQ 1 AND $
                 STRCMP(STRTRIM(trig_type[t],2),'relative') EQ 1 THEN BEGIN 
                 dim=trig_dim[t]+med_mag
              ENDIF 
           ENDIF ELSE dim=med_mag+1.0
           
           ;;bright=0.0
           IF t NE -1 AND STRCMP(STRTRIM(trig_unit[t],2),'mag') EQ 1 THEN BEGIN
              IF STRCMP(STRTRIM(trig_type[t],2),'fixed') EQ 1 THEN BEGIN 
                 bright=trig_bright[t] 
              ENDIF ELSE bright=med_mag-trig_bright[t]
           ENDIF ELSE bright=med_mag-1.0 
           
           ;;--------------------------------------
           ;; Make light curve
           ;;--------------------------------------
    
    plot, [0], [0], /nodata, $
	XR=[MIN(r_jd[r_mytarg])-10.0, MAX(r_jd[r_mytarg])+10.0], $
	YR=[MAX([r_blaz_abs[r_mytarg], r_star_abs[r_mytarg[a]],dim])+0.5,$
    MIN([r_blaz_abs[r_mytarg], r_star_abs[r_mytarg[a]],bright])-0.5], $
	xtitle='MJD', ytitle='Blazar apparent R-Mag.', psym=1, $
	XTICKFORMAT='(I)', TITLE=strtrim(target[tar],2)


    ;;Output colored columns showing previously confirmed triggers
	IF n_prev_trig GT 0 THEN BEGIN 
           FOR u=0,  n_prev_trig-1 DO BEGIN
              OPLOT, [td_mjd[prev_trig[u]], td_mjd[prev_trig[u]]], [-1000,1000], $
                     THICK=15.0, COLOR=FSC_COLOR('Green', /CHECK_CONNECTION)
           ENDFOR
        ENDIF
        ;;Plot any MJD of spectra found
        spec_found = WHERE(specmjd_obj EQ target[tar], n_spec)
        
        IF n_spec GT 0 THEN BEGIN
           FOR u=0, n_spec-1 DO BEGIN
              OPLOT, [specmjd_mjd[spec_found[u]], specmjd_mjd[spec_found[u]]], [-1000, 1000], THICK=15.0,  $
                     COLOR=FSC_COLOR('Slate Blue', /CHECK_CONNECTION), LINE=2
           ENDFOR
        ENDIF
        
        ;;Overplot the reference star
        oploterror, r_jd[r_mytarg[a]], r_star_abs[r_mytarg[a]], r_star_errabs[r_mytarg[a]],$
                    COLOR=FSC_COLOR('Pink', /CHECK_CONNECTION), $
                    ERRCOLOR=FSC_COLOR('Pink', /CHECK_CONNECTION), psym=6
        ;;plot the median and +/-1 mag line	
        oplot, [-100,1D5], [med_mag,med_mag], line=0
	oplot, [-100,1D5], [p1mag,p1mag], line=1
        oplot, [-100,1D5], [m1mag,m1mag], line=1
        ;;Overplot the Blazar
        oploterror, r_jd[r_mytarg], r_blaz_abs[r_mytarg], r_blaz_errabs[r_mytarg], psym=1
        x_psclose
        
        ;;Create claibrated light-curves that report flux density instead of magnitudes
        IF KEYWORD_SET(flux) THEN BEGIN 
           x_psopen, getenv('AUTOM_DIR')+'/plots/'+strtrim(target[tar],2)+'_variab_fluxdensity.ps', /maxs
           ;;Convert blazar magnitudes to flux densities
           flux = MAG_TO_FLUX(r_blaz_abs[r_mytarg]) 
           ;;Convert blazar mag errors to flux densities
           flux_err = MAG_TO_FLUX(r_blaz_abs[r_mytarg] - r_blaz_errabs[r_mytarg]) - flux 
           median_flux = MAG_TO_FLUX(med_mag)
           dim_flux_trigger = MAG_TO_FLUX(trig_dim[t])
           bright_flux_trigger = MAG_TO_FLUX(trig_bright[t])
           ploterror, r_jd[r_mytarg], flux, flux_err, $
                      XR=[MIN(r_jd[r_mytarg])-10.0, MAX(r_jd[r_mytarg])+10.0], $
                      YR=[MIN(flux)*0.9, MAX(flux)*1.1], $
                      xtitle='MJD', ytitle='Flux Density [Jy]', psym=1, $
                      XTICKFORMAT='(I)'
           ;;plot the median and +/-1 mag line	
           plot, [-100,1D5], [median_flux, median_flux], line=0
           oplot, [-100,1D5], [dim_flux_trigger, dim_flux_trigger], line=1
           oplot, [-100,1D5], [bright_flux_trigger, bright_flux_trigger], line=1
           x_psclose
        ENDIF
        ;;Ouput calibrated light-curve into results
        ;;Format is mjd, mag, stddev
	
        ;;Print star calibrated data   
        
        sorted_blazar_jd = sort(r_jd[r_mytarg])
        WRITECOL, getenv('AUTOM_DIR')+'/lightcurve/'+target[tar]+'_calibrated.dat',$
                  r_jd[r_mytarg[sorted_blazar_jd]], $
                  r_blaz_abs[r_mytarg[sorted_blazar_jd]],  r_blaz_errabs[r_mytarg[sorted_blazar_jd]],$
                  FMT='(F7.1, X, F5.2, X, F5.2)'
        
        sorted_star_jd = sort(r_jd[r_mytarg[a]])
        WRITECOL, getenv('AUTOM_DIR')+'/lightcurve/'+target[tar]+'_calibrated_star.dat', $
                  r_jd[r_mytarg[a[sorted_star_jd]]], $
                  r_star_abs[r_mytarg[a[sorted_star_jd]]],  r_star_errabs[r_mytarg[a[sorted_star_jd]]],$
                  FMT='(F7.1, X, F5.2, X, F5.2)'


        ;;***********************************
        ;; PRINT OUT OUTLIERS in lightcurves
        ;;***********************************
        IF KEYWORD_SET(check_outliers) THEN BEGIN
           ;;    SPAWN, "date -d -today +'%y%m%d'", today
        
           ;;in some cases the blazar and star mags are the same...human error somehow
        good_star_median=MEDIAN(r_star_abs[r_mytarg[a[where(ABS(r_star_abs[r_mytarg[a]]-r_blaz_abs[r_mytarg]) $
                                                            GE 0.001)]]])
        ;search for data points that are 3 sigma away from the median values
        ;blaz_outliers=where(ABS((r_blaz_abs[r_mytarg]-med_mag))/r_blaz_errabs[r_mytarg] ge 3.0)
        
        ;in some cases the blazar and star mags are the same...human error somehow
        star_outliers=WHERE(ABS(r_star_abs[r_mytarg[a]]-good_star_median) GE 0.3 $
                            AND ABS(r_star_abs[r_mytarg[a]]-r_blaz_abs[r_mytarg]) GE 0.001)
        
        get_juldate, juldate
        ;;STOP
        response=''
        IF star_outliers[0] NE -1 THEN BEGIN
           READ, response, PROMPT='Open images from '+STRTRIM(target[tar],2)+' Y/n? ' 
           IF STRCMP(response,'Y') EQ 1 OR STRCMP(response,'y') EQ 1 THEN BEGIN
              printf, 6, '#***'+target[tar]+'***'
              printf, 6, '#'+ra[tar]+'    '+dec[tar] 
              printf, 6, 'Blazar outliers'
              printf, 6, '#Star outliers'
              printf, 6, '#Median magnitude='+STRTRIM(good_star_median,2)
              FOR ii=0,N_ELEMENTS(star_outliers)-1 DO BEGIN
                 days_imaged_ago=FLOOR(juldate-2400000.5-r_jd[r_mytarg[a[star_outliers[ii]]]])
                 SPAWN, "date -d "+STRTRIM(days_imaged_ago+1,2)+"-days-ago '+%y%m%d'", yymmdd
                 yymmdd=yymmdd[0]
                 SPAWN, 'ds9 'getenv('AUTOM_DIR')+'/nights/'+STRTRIM(yymmdd,2)+$
                        '/Sci/sci_'+STRTRIM(target[tar],2)+'_R.fits.gz -zoom to fit -zscale'
                 significance_of_outlier=ABS(r_star_abs[r_mytarg[a[star_outliers[ii]]]] - $
                                             MEDIAN(r_star_abs[r_mytarg[a]]))/ $
                                         r_star_errabs[r_mytarg[a[star_outliers[ii]]]]
                 moonpos, r_jd[r_mytarg[a[star_outliers[ii]]]]+2400000.5,moon_ra,moon_dec
                 moon_position=adstring(moon_ra,moon_dec,1)
                 mphase, r_jd[r_mytarg[a[star_outliers[ii]]]]+2400000.5, moon_phase
                 printf, 6, STRTRIM(r_jd[r_mytarg[a[star_outliers[ii]]]],2)+'    '+STRTRIM(yymmdd,2)+'   '+$
                         STRTRIM(r_star_abs[r_mytarg[a[star_outliers[ii]]]],2)+'    '+$
                         STRTRIM(r_star_errabs[r_mytarg[a[star_outliers[ii]]]],2)+'    '+$
                         STRTRIM(significance_of_outlier,2)+'    '+moon_position+'   '+$
                         STRTRIM(moon_phase,2)+'   '+STRTRIM(r_star_rel[r_mytarg[a[star_outliers[ii]]]],2)+'     '+$
                         STRTRIM(r_zp[r_mytarg[a[star_outliers[ii]]]],2)
              ENDFOR
           ENDIF
           printf, 6, ''
           printf, 6, ''
        ENDIF
        
     ENDIF
        
        
        ;;Output relative light-curves
        IF  N_ELEMENTS(r_relative_mytarg) GT 1 THEN BEGIN
           med_mag=median(r_blaz_rel[r_relative_mytarg]-r_star_rel[r_relative_mytarg])
           ;;If there exists a trigger level set by the user
           IF t NE [-1] THEN BEGIN                             
              ;;If trigger type is relative then use +/- from median for tirggers;
              IF trig_type[t] EQ 'relative' AND t NE [-1] THEN BEGIN 
                 p1mag = med_mag + trig_dim[t]
                 m1mag = med_mag - trig_bright[t]
              ENDIF 
              ;;If trigger type is fixed then used the fixed trigger levels,	
              IF trig_type[t] EQ 'fixed' AND t NE [-1] THEN BEGIN 
                 ;;Subtract calibrated median to keep triggers relative
                 m1mag =  med_mag + trig_bright[t] - calibrated_median_mag 
                 p1mag =  med_mag + trig_dim[t] - calibrated_median_mag
              ENDIF
           ENDIF ELSE BEGIN     ;No trigger levels are set so go to defaults of +/- 1 mag from the mean
              p1mag = med_mag + 1.0
              m1mag = med_mag - 1.0
           ENDELSE
           
           x_psopen, getenv('AUTOM_DIR')+'/plots/'+strtrim(target[tar],2)+'_variab.ps', $
                     /maxs      ;Plot the relative light-curve
           ;;Add errors in quadrature
           relative_error = SQRT(r_blaz_errrel[r_relative_mytarg]^2 + r_star_errrel[r_relative_mytarg]^2) 
           
           PLOT, [0], [0], /NODATA, xtitle='MJD', $
                 ytitle=Textoidl('m_{blazar}-m_{star}'), psym=1, $
                 xrange=[MIN(r_jd[r_relative_mytarg])-10., MAX(r_jd[r_relative_mytarg])+10], $
                 yrange=[MAX(r_blaz_rel[r_relative_mytarg]-r_star_rel[r_relative_mytarg])+0.5,$
                         MIN(r_blaz_rel[r_relative_mytarg]-r_star_rel[r_relative_mytarg])-0.5],$
                 XTICKFORMAT='(I)'
           IF n_prev_trig GT 0 THEN BEGIN ;Output colored columns showing previously confirmed triggers
              FOR u=0,  n_prev_trig-1 DO BEGIN
                 OPLOT, [td_mjd[prev_trig[u]], td_mjd[prev_trig[u]]], [-1000,1000], THICK=15.0, $
                        COLOR=FSC_COLOR('Green', /CHECK_CONNECTION)
              ENDFOR
           ENDIF
           ;;Plot any MJD of spectra found
           spec_found = WHERE(specmjd_obj EQ target[tar], n_spec)
           IF n_spec GT 0 THEN BEGIN
              FOR u=0, n_spec-1 DO BEGIN
                 OPLOT, [specmjd_mjd[spec_found[u]], specmjd_mjd[spec_found[u]]], [-1000, 1000], THICK=15.0,  $
                        COLOR=FSC_COLOR('Slate Blue', /CHECK_CONNECTION), LINE=2
              ENDFOR
           ENDIF
           oploterror, r_jd[r_relative_mytarg], r_blaz_rel[r_relative_mytarg]-r_star_rel[r_relative_mytarg],$
                       relative_error, psym=1
           ;;plot the median and +/-1 mag 
           oplot, [-100,1D5], [med_mag,med_mag], line=0
           oplot, [-100,1D5], [p1mag,p1mag], line=1
           oplot, [-100,1D5], [m1mag,m1mag], line=1
           
	x_psclose
        
        ;;Output zero piont light-curve
        x_psopen, getenv('AUTOM_DIR')+'/plots/'+strtrim(target[tar],2)+'_variab_zeropoint.ps', $
                  /maxs         ;Plot the relative light-curve
        PLOT, [0], [0], /NODATA, xtitle='MJD', $
              ytitle='Zero Point', psym=1,$
              xrange=[MIN(r_jd[r_mytarg])-10., MAX(r_jd[r_mytarg])+10], $
              yrange=[MIN(r_zp[r_mytarg])-0.5, MAX(r_zp[r_mytarg])+0.5],$
              XTICKFORMAT='(I)'
        IF n_prev_trig GT 0 THEN BEGIN ;Output colored columns showing previously confirmed triggers
           FOR u=0,  n_prev_trig-1 DO BEGIN
              OPLOT, [td_mjd[prev_trig[u]], td_mjd[prev_trig[u]]], [-1000,1000], $
                     THICK=15.0, COLOR=FSC_COLOR('Green', /CHECK_CONNECTION)
           ENDFOR
        ENDIF
        ;;Plot any MJD of spectra found
        spec_found = WHERE(specmjd_obj EQ target[tar], n_spec)
        IF n_spec GT 0 THEN BEGIN
           FOR u=0, n_spec-1 DO BEGIN
              OPLOT, [specmjd_mjd[spec_found[u]], specmjd_mjd[spec_found[u]]], [-1000, 1000], THICK=15.0,  $
                     COLOR=FSC_COLOR('Slate Blue', /CHECK_CONNECTION), LINE=2
           ENDFOR
        ENDIF
        oploterror, r_jd[r_mytarg], r_zp[r_mytarg], r_zp_err[r_mytarg], psym=1
        x_psclose	
        
     endif                      

        ;;Ouput relative light-curve into results
        ;;Format is mjd, mag, stddev
	WRITECOL, getenv('AUTOM_DIR')+'/data/results/'+target[tar]+'_relative.dat', $
                  r_jd[r_relative_mytarg],  r_blaz_rel[r_relative_mytarg]-r_star_rel[r_relative_mytarg],$
                  relative_error, FMT='(F7.1, X,F5.2, X,F5.2)'
     endif
     endif
     
     ;;Update test light-curve webpage
     blz_name = strtrim(target[tar],2)
     
     ned_blz_name = repstr(blz_name,'+','%2B')
     ned_blz_name = repchr(ned_blz_name,'_','+')
     
     printf, 12, '<a href="lightcurve_webpage/'+blz_name+'.html">'+repchr(blz_name,'_',' ')+'</a><br><br>'
     printf, 5, '<th>'+repchr(blz_name,'_',' ')+'</th>'
     printf, 5, '<tr><td>'
     IF n_elements(r_mytarg) GT 1 THEN BEGIN ;If light-curves exist
        printf, 5, '<a href="images/'+blz_name+'_variab_calib.jpg">'
        printf, 5, '<img src="images/'+blz_name+'_variab_calib.jpg" width="450" height="350" border="0"></a></td>'
     ENDIF ELSE printf, 5, 'Not enough data to make light-curve.</td>' ;No light-curve found
     printf, 5, '<td>Name: '+repchr(blz_name,'_',' ')+'<br>'
     printf, 5, '<font size="-1">'
     printf, 5, 'Further Information:'
     printf, 5, '<a href="http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?objname='+ned_blz_name+$
             '&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox='+$
             'J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES">NED</a><br>'
     printf, 5, '</font>'
     printf, 5, '<br> RA: '+ra[tar]+'<br> Dec: '+dec[tar]+'<br><br>'
     printf, 5, '<font size="-1">'
     printf, 5, 'Observable by <a href="http://tevcat.uchicago.edu/visplot.cgi?tname='+blz_name+' &ra='+ra[tar]+' &dec='+dec[tar]+'&lat=31.68 &date='
     printf, 5, ut_day+'-'+ut_month+'-'+ut_year+' &lon=-110.86&mode=1">VERITAS</a>, '
     printf, 5, '<a href="http://tevcat.uchicago.edu/visplot.cgi?tname='+blz_name+' &ra='+ra[tar]+' &dec='+dec[tar]+'&lat=19.0 &date='
     printf, 5, ut_day+'-'+ut_month+'-'+ut_year+' &lon=-155.47 &mode=1">Keck</a>, '
     printf, 5, '<a href="http://tevcat.uchicago.edu/visplot.cgi?tname='+blz_name+' &ra='+ra[tar]+' &dec='+dec[tar]+'&lat=37.34311 &date='
     printf, 5, ut_day+'-'+ut_month+'-'+ut_year+' &lon=-121.64278 &mode=1">Lick</a>?<br>'
     printf, 5, '<font size="-2">(generated by <a href="http://tevcat.uchicago.edu/CustomVis.pl">TeVCat</a>)</font><br><br>'
     printf, 5, '</font>'
     days_ago = CURRENT_MJD() - current_jd
     printf, 5, 'Latest imake taken on MJD:<br>'+strmid(strcompress(current_jd, /remove_all),0,7)
     printf, 5, 'which was '+STRING(days_ago, FORMAT='(F10.1)')+' days ago.<br><br>'
     
     med_mag = LIGHT_CURVE_MEDIAN_MAG(blz_name)
                                ;print median magnitudes and trigger levels to webpage and write to file med_mag.txt
     IF (med_mag NE -1) THEN BEGIN
        IF trig_type[t] EQ 'relative' THEN BEGIN ;If trigger type is relative then use +/- from median for tirggers
           highstatetrig=med_mag-trig_bright[t]
           lowstatetrig=med_mag+trig_dim[t]
        ENDIF ELSE BEGIN        ;If trigger type is fixed then used the fixed trigger levels, ignoring the median
           lowstatetrig=trig_dim[t] 
           highstatetrig=trig_bright[t]
        ENDELSE
        printf, 5, 'High state threshold: '+STRTRIM(highstatetrig,2)+'<br>'
        printf, 5, 'Low state threshold: '+STRTRIM(lowstatetrig,2)+'<br><br>'
        printf, 5, 'Median Magnitude: '+STRING(med_mag)+' <br>'
        printf, 10, FORMAT='(A, 5X, F0.3, 5X, F0.3, 5X, F0.3)', blz_name, med_mag, lowstatetrig, highstatetrig
        last_point = r_blaz_abs[r_mytarg[N_ELEMENTS(r_mytarg)-1]]
        printf, 5, 'Newest deviation from median: '+STRING(ABS(med_mag-last_point))+' <br><br>'
     ENDIF ELSE BEGIN
        printf, 5, 'Median Magnitude: Unavailable <br>'
        printf, 10, FORMAT='(A, 5X, A, 5X, A, 5X, A)', blz_name, 'Unavailable', 'Unavailable', 'Unavailable'
        printf, 5, 'Newest deviation from median: Unavailable <br><br>'
     ENDELSE
     
     ;;print subclass and log nu_syn to webpage
     index = subcl_and_lognu(blz_name,names_array)
     IF (index NE -1) THEN BEGIN
        printf, 5, 'Subclass: '+subclass[index]+' <br>'
        IF(lognu[index] EQ -1) THEN BEGIN
           printf, 5, 'Log nu_syn: Unknown <br>'
        ENDIF ELSE BEGIN 
           printf, 5, 'Log nu_syn: '+lognu[index]+' <br><br>'
        ENDELSE
     ENDIF ELSE BEGIN
        printf, 5, 'Subclass: Unavailable <br>'
        printf, 5, 'Log nu_syn: Unavailable <br><br>'
     ENDELSE
     
     printf, 5, '<a href="images/'+blz_name+'_variab.jpg">Relative Light-Curve</a><br><br>'
     printf, 5, '<a href="images/'+blz_name+'_variab_zeropoint.jpg">Zero Point Light-Curve</a><br><br>'
     printf, 5, '<a href="Finder_Charts/'+blz_name+'_fndchr.jpg">Finder Chart</a><br><br>'
     printf, 5, '<a href="trigger_webpage/thumbnails/'+blz_name+'.jpg">Latest Super-LOTIS image</a><br><br>'
     printf, 5, '<a href="results_calibrated/'+blz_name+'_calibrated.dat">Calibrated photometry data</a><br><br>'
     printf, 5, '<a href="results_relative/'+blz_name+'_relative.dat">Relative photometry data</a><br><br>'
     printf, 5, '</td></tr>'
     
     ;;creates light curve page for each blazar
     
     automated_mklcpage, r_mytarg, ut_year, ut_month, ut_day, blz_name, ned_blz_name,$
                           strcompress(current_jd, /remove_all), ra[tar], dec[tar], days_ago,$
                           med_mag, highstatetrig, lowstatetrig, last_point, index, subclass,$
                           lognu, significance, alert 
     ;;reset alert for next blazar
    alert=0
 endfor
  
  

  ;;Close all the output files
  CLOSE, 3                      ;Close text file holding trigger results
  CLOSE, 4                      ;Close webpage for triggers
  printf, 5, '</table>'
  CLOSE, 5                      ;Close test webpage for light-curves
  IF keyword_set(check_outliers) THEN CLOSE, 6
  PRINTF, 12, '</body></html>'
  CLOSE, 12
  
  ;;copies the median mag and trigger levels into each day
  dateString='date -d -1-days +'
  dateString=dateString+"'%y%m%d'"
  SPAWN, dateString, dateOutput
  dateOutput=STRTRIM(STRING(dateOutput[0]),2)
  CLOSE, 10                     ;Close median_mag.txt
  SPAWN, 'cp  '+getenv('AUTOM_DIR')+'/data/median_mag.txt '+getenv('AUTOM_DIR')+'/data/'+dateOutput+$
         '/'+dateOutput+'_trigger_levels.txt'
  
  ;;Update webpage
  IF KEYWORD_SET(update_webpage) THEN BEGIN
     SPAWN, 'scp '+getenv('AUTOM_DIR')+'/data/'+dateOutput+'/'+dateOutput+$
            '_trigger_levels.txt kkaplan@vhe3.ucsc.edu:/home/vhep/afurniss/public_html/private/BlazarMonitoring/trigger_levels'
     CD, 'PLOTS'
     convert_pstojpg

     ;;;;;;;; ---- THIS PART WILL BREAK AND NEEDS TO BE UPDATED  ---- ;;;;;

     SPAWN, 'scp *variab*jpg kkaplan@vhe3.ucsc.edu:/home/vhep/afurniss/public_html/private/BlazarMonitoring/images'
                                ;Update test trigger webpage automatically
     SPAWN, 'scp -r /b/Blazars/data/trigger_webpage/* kkaplan@vhe3.ucsc.edu:/home/vhep/afurniss/public_html/private/BlazarMonitoring/trigger_webpage'
     SPAWN, 'scp -r /b/Blazars/data/lightcurve_webpage/* kkaplan@vhe3.ucsc.edu:/home/vhep/afurniss/public_html/private/BlazarMonitoring/lightcurve_webpage'
     SPAWN, 'scp  /b/Blazars/data/Light-Curves.html kkaplan@vhe3.ucsc.edu:/home/vhep/afurniss/public_html/private/BlazarMonitoring/'
     SPAWN, 'scp /b/Blazars/data/results/*calibrated.dat kkaplan@vhe3.ucsc.edu:/home/vhep/afurniss/public_html/private/BlazarMonitoring/results_calibrated/'
     SPAWN, 'scp /b/Blazars/data/results/*relative.dat kkaplan@vhe3.ucsc.edu:/home/vhep/afurniss/public_html/private/BlazarMonitoring/results_relative/'	
     SPAWN, 'scp /b/Blazars/data/recent_observations.html kkaplan@vhe3.ucsc.edu:/home/vhep/afurniss/public_html/private/BlazarMonitoring/recent_observations/'
     SPAWN, 'scp /b/Blazars/data/All-Light-Curves.html kkaplan@vhe3.ucsc.edu:/home/vhep/afurniss/public_html/private/BlazarMonitoring/'
     CD,'..'
     thumbnail_page
  ENDIF
  
   
  ;;Make a plot showing the times a blazar is imaged
  X_RADEC, ra, dec, deg_ra, deg_dec ;Convert RA and Dec. from sexigasimal to decimal degrees
  x_PSOPEN, getenv('AUTOM_DIR')+'/plots/slotis_image_timing.ps', /MAXS
  PLOT, [0], [0], /NODATA, XR=[55190, MAX(r_jd)+10], YR=[0,360], XTITLE='MJD', YTITLE='RA (deg.)', XTICKFORMAT='(I)'
  FOR i=0, ntarget-1 DO BEGIN
     j = WHERE(r_targ EQ target[i] and r_jd NE 0) ;Find all instances of a target in the photometry data
     IF j NE [-1] THEN BEGIN
        ra_array = MAKE_ARRAY(N_ELEMENTS(j), VALUE=deg_ra[i]) ;Create an array for RA that will be the y for the plot
        OPLOT, r_jd[j], ra_array, PSYM=1
     ENDIF
  ENDFOR
  x_psclose
  
  
end
