;+
; NAME:
;   automated_work
;
; PURPOSE:
; Procedure that runs the superlotis pipeline and produce light curves for object
;
; 1) Determine if a night has already been reduced and analyzed
; 2) If not, then retrieve the data
; 3) Run data reduction and analysis on new data
; 4) Clean files storing only raw and fully reduced, if desired
;
;    
; CALLING SEQUENCE:
;
; automated_work, date, REQUEST=, /CLEAN,/NOGETDATA, $
;       /BADSEEING, /CLOBBER
;
; INPUTS:
;   
; date        -- string day to process  
; request     -- the queue file. File format
;                
; 
; OPTIONAL INPUTS:
;  
; /clean       -- if set deletes intermediate step files 
; /nogetdata   -- do not download data again
; /badseeing   -- increase the aperture size used for photometry
; /clobber     -- re-reduced the data overwriting everything (as in
;                it is clobbering time).
;
; OUTPUTS: 
;
; Creates a bunch of reduced science frame  
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; Only one filter is read with this current version (R band).
;
; EXAMPLES:
;
; Run the night 04 July 2012 which has not been retrived and 
; clean the directory of intermediate step files.
;       automated_work, '120704', REQUEST='path/to/file', /CLEAN
;
; Rerun 04 July 2012 which the data has already been retrieved 
; and clean the directory of intermediate step files.
;       automated_work, '120704', REQUEST='path/to/file', /CLEAN, $
;       /NOGETDATA, /CLOBBER
;
; BUGS:
;
; PROCEDURES CALLED:
;   automated_getdata, automated_checklognames, automated_cleanlog, 
;   automated_ccdproc, automated_makeplan, automated_stackimages,
;   automated_phot
;
; REVISION HISTORY:
;    05-Jun-2012  Written and revised by MF
;
;-
;------------------------------------------------------------------------------
;; 


pro automated_work, date, request=request, clean=clean, nogetdata=nogetdata, $
                    badseeing=badseeing, clobber=clobber, FLWO=FLWO

   ;;journal, date+'_reduction_log'
   splog, 'Starting at ', SYSTIME()
  
   ;;-----------------------------------------
   ;;GET DATA 
   ;;------------------------------------------  

   ;;If data folder does not exist, make a directory for the night
   ;;otherwise remove extra files and re-reduce
   
   

   spawn, 'mkdir -p nights/'+date ;create data folder for this night
   spawn, 'cp '+request+' '+date ;Copy target list into data folder
   if (keyword_set(clobber) && ~KEYWORD_SET(FLWO)) then begin
     ; spawn, 'rm nights/'+date+'/Sci/*' ;Remove science files
     ; spawn, 'rm nights/'+date+'/*fit*' ;Remove fits files
     ; spawn, 'rm nights/'+date+'/*.str*' ;Remove stray log files
      spawn, 'tar -zxvf nights/'+date+'/slotis*'+' -C nights/'+date ;unzip the .tar file
      spawn, 'gunzip nights/'+date+'/*fit*' ;unzip the data files
   endif


   if (keyword_set(clobber) && KEYWORD_SET(FLWO)) then begin
     ; spawn, 'rm nights/'+date+'/Sci/*' ;Remove science files
     ; spawn, 'rm nights/'+date+'/*fit*' ;Remove fits files
     ; spawn, 'rm nights/'+date+'/*.str*' ;Remove stray log files
      spawn, 'tar -zxvf nights/'+date+'/FLWO_1*'+' -C nights/'+date ;unzip the .tar file
      spawn, 'gunzip nights/'+date+'/*fit*' ;unzip the data files
      
   endif




   ;;Get into data foldersl
   splog, 'Entering ./nights/'+date
   cd, './nights/'+date, current=main_dir 
  
   ;;get the data if needed
   if ~keyword_set(nogetdata) and ~keyword_set(clobber) then begin
      ;;check if files are already extracted
      file=file_info('bias1.fits')
      if(file.exists eq 1) then begin 
         splog, 'Data file already extracted' 
      endif else begin
         automated_getdata, date
      endelse
   endif
  

   ;;-----------------------------------------
   ;;REDUCTION 
   ;;------------------------------------------  

   ;;check log and replace old names with new names
   IF ~KEYWORD_SET(FLWO) THEN BEGIN

      automated_checklognames, 'obs_'+date+'.log' 
      splog, 'Starting Slotis reduction...'   
      ;;Clean logfile
      automated_cleanlog, date, request
   ENDIF ELSE BEGIN
      splog, 'Starting FLWO reduction...'
      automated_makelog,date

   ENDELSE

   ;;Makeplan
   IF ~KEYWORD_SET(FLWO) THEN BEGIN
      automated_makeplan, 'obs_'+date+'.log'
   ENDIF ELSE BEGIN
      automated_makestr, 'FLWO_obs_'+date+'.log'
   ENDELSE
  
   ;;Call ccdproc
   IF ~KEYWORD_SET(FLWO) THEN BEGIN   
       automated_ccdproc, 'obs_'+date+'.log.str'
   ENDIF ELSE BEGIN
       logstr = STRARR(4)
       logstr[0] = 'FLWO_obs_'+date+'.log_1.str'
       logstr[1] = 'FLWO_obs_'+date+'.log_2.str'
       logstr[2] = 'FLWO_obs_'+date+'.log_3.str'
       logstr[3] = 'FLWO_obs_'+date+'.log_4.str'
       automated_FLWOccdproc,logstr
   ENDELSE
  
   ;;Make final frames
 ;  IF ~KEYWORD_SET(FLWO) THEN BEGIN
 ;  automated_stackimages, 'obs_'+date+'.log.str', /wcs, badseeing=badseeing
 ;  ENDIF ELSE BEGIN
      
      automated_stackimages,logstr[0],/wcs,badseeing=badseeing
      
 ;  ENDELSE
  
   ;;Enter science directory
   splog, 'Entering ./Sci/'
   cd, './Sci/'
  
   ;;Autophot
   filter=['r','VH','i','BH']
   FOR i=0,3 DO BEGIN
      automated_phot, date, FLWO=FLWO,filter[i] ;;Run with R band for now. Could expand this 
   ENDFOR

   ;;--------------------------------------  
   ;;CLEAN 
   ;;--------------------------------------
        
;   splog, 'Compress...(This take a while, but then you are done!)'
;  
;   if keyword_set(clean) then begin
;      splog, 'Clean some stuff'
;      ;;Clean stuff
;      splog, 'Entering ./'+date
;      cd, './'+date 
;      spawn, 'rm -f *.fits*'
;      spawn, 'rm -f *.fit*'
;      spawn, 'rm -f obs_'+date+'.log.str*'
;      spawn, 'gzip *.tar'
;      ;;Compress science images
;      splog, 'Entering ./Sci/'
;      cd, './Sci/' 
;      spawn, 'gzip sci_*.fits'
;      ;;Back to main dir
;      splog, 'Entering '+main_dir
;      cd, main_dir
;   endif else begin
;      ;;Compress reduction
;      splog, 'Entering ./'+date
;      cd, './nights/'+date 
;      spawn, 'gzip *.fits'
;      spawn, 'gzip *.fit'
;      spawn, 'gzip *.tar'
;      ;;Compress science images
;      splog, 'Entering ./Sci/'
;      cd, './Sci/' 
;     
;      ;;this was added by Matt. I am commenting this out for now, as  
;      ;;I think there are some mismatched keywords. Also, the code
;      ;;should almost always run with /clean keyword set.
;      ;;thumbnail_page, date, /already_made  
;      spawn, 'gzip sci_*.fits'
;      ;;back to main dir
;      splog, 'Entering '+main_dir
;      cd, main_dir
;   endelse
;  
   splog, 'All done at ', SYSTIME()
   journal
  
end
