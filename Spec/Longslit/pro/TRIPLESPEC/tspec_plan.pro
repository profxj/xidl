;+
; NAME:
;   tspec_plan
;
; PURPOSE:
;   Create plan file(s) for running the GNIRS pipeline
;
; CALLING SEQUENCE:
;   tspec_plan, [ indir, fileexpr, planfile= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   fileexpr   - File names in the input directory; default to '*.fits*'
;   indir      - Input directory(s) for reading files;
;                default to current directory
;   planfile   - Output plan file; default to 'plan.par'.
;                This file is put in the same directory as the raw data files.
;
; OUTPUT:
;
; COMMENTS:
;   One plan file is made for each input directory.
;
;   The following flavors of images are listed:
;     ACQ
;     FLAT
;     ARC
;     TELL
;     SCIENCE
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   fileandpath()
;   headfits()
;   idlutils_version()
;   splog
;   sxpar()
;   yanny_write
;
; INTERNAL SUPPORT ROUTINES:
;   tspecs_plan_struct()
;
; REVISION HISTORY:
;   Dec-2012  Written by JXP; Modified gnris_plan by JFH
;-
;------------------------------------------------------------------------------
FUNCTION tspec_plan_struct, nfile
planstr = create_struct(name = 'lexp' $
                        , 'FILENAME', ''  $  
                        , 'GEMINDX', 0L   $  
                        , 'FLAVOR', ''    $
                        , 'TARGET', ''    $
                        , 'EXPTIME', 0.0  $
                        , 'GRATWAVE', 0.0 $
                        , 'GRATTILT', 0.0 $
                        , 'GRATORD', 0L   $
                        , 'SLIT',    ''   $
                        , 'GROUP', 0L   $
                        , 'AIRMASS', 0.0  $
                        , 'SEQ_NUM', 0L   $
                        , 'ABBA_NUM', 0L)

return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro tspec_plan, indir, fileexpr, planfile = planfile, tell_time = tell_time, $
                USER_DEFINE=user_define

   COMMON SITE, lat, lng, tzone
   DRADEG = 180.d0/!dpi
   DAY_SEC = 24.0D*60.0D*60.0D
   
   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(fileexpr)) then fileexpr = '*.fits*'
   if (NOT keyword_set(planfile)) then planfile = 'plan.par'
   IF NOT KEYWORD_SET(TELL_TIME) THEN TELL_TIME = 10.0

   if not keyword_set(INDIR) then begin
       spawn, '\ls -d '+indir, dirlist
       if (NOT keyword_set(dirlist)) then begin
           splog, 'No input directories found'
           return
       endif
       ndir = n_elements(dirlist)
   endif else begin
       dirlist = [indir]
       ndir = n_elements(dirlist)
   endelse
   
   ;----------
   ; Loop over each input directory

   ;; USER-DEFINE
   if keyword_set(USER_DEFINE) then begin
      readcol, USER_DEFINE, user_file, user_type, format='A,A'
   endif

   for idir = 0L, ndir-1L do begin
       splog, 'Working on directory ', dirlist[idir]
       cd, dirlist[idir], current = olddir
       if (idir EQ 0) then origdir = olddir
       filenames = findfile(fileexpr, count = nfile)
       splog, 'Number of FITS files found: ', nfile
       
       if (nfile GT 0) then begin
           planstr = tspec_plan_struct(nfile)
           planstr.filename = filenames
           mjd_obs = dblarr(n_elements(planstr))
           date_obs = strarr(n_elements(planstr))
           for i = 0L, nfile-1L do begin
               hdr = headfits(filenames[i])
               observ = strtrim(sxpar(hdr, 'OBSERVAT'),2) ;; APO vs Palomar
               IF i EQ 0 THEN mjd_ref = double(sxpar(hdr, 'MJD_OBS'))
               if (size(hdr, /tname) EQ 'STRING') then begin
                   acqmir =  strtrim(sxpar(hdr, 'ACQMIR'))
                   slit   =  strcompress(strtrim(sxpar(hdr, 'SLIT'),2),/remove_all)
                   exptime = float(strtrim(sxpar(hdr, 'EXPTIME')))
                   airmass = float(strtrim(sxpar(hdr, 'AIRMASS')))
                   ;planstr[i].GEMINDX =  gemindx_extract(filenames[i])
                   planstr[i].EXPTIME = exptime

                   case observ of
                      'APO': begin
                         planstr[i].slit = slit
                         planstr[i].TARGET = strtrim(sxpar(hdr, 'OBJNAME'))
                         obstype = strtrim(sxpar(hdr, 'IMAGETYP'))
                         mjd = strmid(sxpar(hdr,'DATE-OBS'), 0, 10)
                         time = strmid(sxpar(hdr,'DATE-OBS'), 12)
                      end
                      else: begin ;; Default to Palomar
                         planstr[i].TARGET = strtrim(sxpar(hdr, 'CATNAME'))
                         obstype = strtrim(sxpar(hdr, 'OBSTYPE'))
                         mjd = strmid(sxpar(hdr,'UTSHUT'), 0, 9)
                         time = sxpar(hdr, 'TIME')
                      end
                   endcase
                   ;planstr[i].TARGET = strtrim(sxpar(hdr, 'OBJECT'))

                   ;planstr[i].SLIT = strtrim(sxpar(hdr, 'SLIT'))
                   ;planstr[i].GRATWAVE = float(sxpar(hdr, 'GRATWAVE'))
                   ;planstr[i].GRATORD  = long(sxpar(hdr, 'GRATORD'))
                   ;planstr[i].GRATTILT = float(sxpar(hdr, 'GRATTILT'))
                   planstr[i].AIRMASS  = airmass
                   ;; mjd  - mjd_reference (first object in list)
                   mjd_obs[i] = x_setjdate(mjd, time)
                   ;date_obs[i]         = sxpar(hdr, 'DATE-OBS')
                   planstr[i].FLAVOR  = obstype 
                   if keyword_set(USER_DEFINE) then begin
                      mt = where(strmatch(user_file, filenames[i]), nmt)
                      case nmt of
                         0: 
                         1: begin
                            planstr[i].flavor = user_type[mt]
                            print, 'tspec_plan: Over riding FLAVOR for ', filenames[i]
                         end
                         else: stop
                      endcase
                   endif
                   ;IF obstype EQ 'FLAT' THEN  planstr[i].FLAVOR = 'FLAT' $
                   ;ELSE IF obstype EQ 'ARC'  THEN  planstr[i].FLAVOR = 'ARC'  $
                   ;ELSE IF (acqmir EQ 'In' OR slit EQ 'Acquisition') $
                   ;  THEN planstr[i].FLAVOR = 'ACQ' $
                   ;ELSE IF obstype EQ 'OBJECT' AND  exptime LE TELL_TIME THEN $
                   ;  planstr[i].FLAVOR = 'TELL' $
                   ;ELSE IF obstype EQ 'OBJECT' AND exptime GE TELL_TIME $
                   ;  THEN planstr[i].FLAVOR = 'SCIENCE' $
                   ;ELSE BEGIN
                   ;    planstr[i].FLAVOR = 'UNKNOWN'
                   ;    splog, 'Unrecognized file'
                   ;ENDELSE
               ENDIF
            ENDFOR

           ;; Zero out BAD files
           bad = where(strmatch(strtrim(planstr.FLAVOR,2), 'BAD'), complement=good, nbad)
           if nbad GT 0 then begin
              print, 'tspec_plan:  Removing bad files!'
              planstr = planstr[good]
              mjd_obs = mjd_obs[good]
              date_obs = date_obs[good]
           endif

;          Identify different groups and label them 0 to ngroup-1
;          at the moment, grouping is done by date. One could also do 
;          have different setups and different groups, but chances
;          are observations on the same date use the same setup for my data
;           setup_list = string(planstr.GRATWAVE, FORMAT = '(F6.4)') + '-' $
;             + string(planstr.GRATORD, FORMAT = '(F6.4)') + '-' $
;             + string(planstr.GRATORD, FORMAT = '(I1.1)') 
           date_uniq = date_obs[uniq(date_obs, sort(date_obs))]
           ngroup = n_elements(date_uniq)
           FOR j = 0L, ngroup-1L DO BEGIN
               this_group = WHERE(date_obs EQ date_uniq[j], nthis)
               planstr[this_group].GROUP = j
           ENDFOR
;          For each setup, loop over the science frames grouping them into
;          ABBA sequences


           FOR j = 0L, ngroup-1L DO BEGIN
               group_inds = WHERE(planstr.GROUP EQ j $
                                  AND (planstr.FLAVOR EQ 'object' OR  $
                                       planstr.FLAVOR EQ 'telluric'), nsci)
               IF nsci GT 0 THEN targets = planstr[group_inds].TARGET $
               ELSE CONTINUE
               uniq_targets = targets[uniq(targets, sort(targets))]
               ntarg = n_elements(uniq_targets)
               For k = 0L, ntarg-1L DO BEGIN
                   this_targ = WHERE(planstr.TARGET EQ uniq_targets[k] $
                                     AND planstr.GROUP EQ j $ 
                                     AND (planstr.FLAVOR EQ 'object' OR  $
                                          planstr.FLAVOR EQ 'telluric'), nabba)
;                  Sort these images with respect to observation time
                   this_mjd = mjd_obs[this_targ]
                   mjd_sort = sort(this_mjd)
                   this_targ = this_targ[mjd_sort]
                   this_mjd = this_mjd[mjd_sort]
                   IF nabba MOD 4 NE 0 THEN splog, $
                     'WARNING: The number of exposures for target:' + $
                     uniq_targets[k] + ' is not divisible by 4'
                   abba_ctr = 0L
                   seq_ctr  = 0L
                   mjd_last = this_mjd[0]
                   FOR ii = 0L, nabba-1L DO BEGIN
                       bigindex = this_targ[ii]
                       IF (ABBA_CTR EQ 0) THEN BEGIN
                           planstr[bigindex].abba_num = abba_ctr
                           planstr[bigindex].seq_num = seq_ctr
                           abba_ctr = abba_ctr+1L
                       ENDIF ELSE IF ABBA_CTR GT 0 THEN BEGIN
                           mjd_diff = this_mjd[ii]-mjd_last
;                          Is this the next step in the same sequence?
                           IF mjd_diff LE (3*planstr[bigindex].EXPTIME/DAY_SEC > 0.002) $
                           THEN BEGIN
;                          If yes then assign the sequence and increment counter
                               planstr[bigindex].abba_num = abba_ctr
                               planstr[bigindex].seq_num = seq_ctr
                               IF (abba_ctr EQ 3) then begin
                                   abba_ctr = 0
                                   seq_ctr  = seq_ctr + 1
                               ENDIF ELSE abba_ctr = abba_ctr + 1
                           ENDIF ELSE BEGIN
                              stop
                               splog, 'Found file out of sequence: ' + $
                                      planstr[bigindex].FILENAME  + $
                                      '  Target: ' + planstr[bigindex].TARGET
                               splog, 'Resetting sequence'
                               abba_ctr = 0L
                               seq_ctr = seq_ctr+1
                               planstr[bigindex].abba_num = abba_ctr
                               planstr[bigindex].seq_num = seq_ctr
                           ENDELSE
                       ENDIF
                       mjd_last = this_mjd[ii]
                   ENDFOR
               ENDFOR
           ENDFOR
           
           logfile = repstr(planfile, '.par', '') + '.log'
           plotfile = repstr(planfile, '.par', '') + '.ps'

           
           hdr = ''
           hdr = [hdr, '# Sequence grouped by GRATWAVE+GRATTILT+TARGET']
           hdr = [hdr, ' ']
           hdr = [hdr, "logfile '" + logfile + "'   # Log file"]
           hdr = [hdr, "plotfile '" + plotfile + "'   # Plot file"]
           hdr = [hdr, "indir '" + dirlist[idir] + "'   # Raw data directory"]
           hdr = [hdr, "tempdir Temp     # Temporary (working) directory"]
           hdr = [hdr, "scidir  Science  # Science output directory"]
           hdr = [hdr, "idlutilsVersion '" + idlutils_version() $
          + "'  # Version of idlutils when building plan file"]
;           hdr = [hdr, "LongslitVersion '" + longslit_version() $
;          + "'  # Version of Longslit when building plan file"]

;        Write this plan file
           cd, olddir
           yanny_write, planfile, ptr_new(planstr), hdr = hdr, /align
       endif
       
   endfor                       
; End loop over directories
   
   cd, origdir

   ;; Make some sub-directories
   QApath = 'QA/'
   IF file_test(qapath, /dir) EQ 0 THEN spawn, 'mkdir ' + qapath
   newpath = 'Flux/'
   IF file_test(newpath, /dir) EQ 0 THEN spawn, 'mkdir ' + newpath
   newpath = 'FSpec/'
   IF file_test(newpath, /dir) EQ 0 THEN spawn, 'mkdir ' + newpath
   
   return
end
;------------------------------------------------------------------------------
