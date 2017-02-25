;+
; NAME:
;   sinfoni_plan
;
; PURPOSE:
;   Create plan file(s) for running the SINFONI pipeline
;
; CALLING SEQUENCE:
;   sinfoni_plan, [ fileexpr, indir, planfile= ]
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
;   niri_plan_struct()
;
; REVISION HISTORY:
;   07-Jun-2006  Written by David Schlegel, LBL.
;-
;------------------------------------------------------------------------------

FUNCTION sinfoni_frames, inds
skystr = '[' + strcompress(string(inds[0]), /rem)
FOR j = 1L, n_elements(inds)-1L DO BEGIN
    skystr = skystr + ', ' + strcompress(string(inds[j]), /rem)
ENDFOR    
skystr = skystr + ']'

RETURN, skystr
END


FUNCTION sinfoni_sequence, seq_num, nseq
  CASE nseq OF
     1: BEGIN 
        ;; calibration image
        sky_ind = [1] 
     END
     2: BEGIN 
        ;; 2-pt Dither pattern, each uses the other
        ;;         1, 2
        sky_ind = [2, 1] 
     END
     4: BEGIN 
        ;; 4-pt Dither pattern
        ;;  none uses 1 image as sky. 3 and 4 use each other
        ;;         1, 2, 3, 4 
        sky_ind = [2, 3, 4, 3] 
     END
     6: BEGIN 
        ;; 6-pt Dither pattern
        ;; 
        ;;  none use 1 image as sky. 5 and 6 use each other
        ;;         1, 2, 3, 4, 5, 6
        sky_ind = [2, 3, 4, 5, 6, 5]
     END  
     8: BEGIN 
        ;;  none uses 0 image as sky. 2 and 3 use each other
        ;;         1, 2, 3, 4, 5, 6, 7, 8
        sky_ind = [2, 3, 4, 5, 6, 7, 8, 7]
     END
     24: BEGIN 
        ;;  none uses 0 image as sky. 2 and 3 use each other
        sky_ind = lindgen(nseq) + 2L
        sky_ind[nseq-1L] = nseq-1
     END  
     ELSE: message, 'Unsupported sequence number. Add it here'
  ENDCASE
  RETURN, sky_ind[seq_num-1]
END


FUNCTION sinfoni_plan_struct, nfile
planstr = create_struct(name = 'lexp' $
                        , 'FILENAME', ''  $  
                        , 'GROUP', -1L     $
                        , 'ISEQ', 1L   $    
                        , 'SEQ_NUM', 0L   $
                        , 'NSEQ', 0L   $  
                        , 'FLAVOR', ''    $
                        , 'TARGET', ''    $
                        , 'EXPTIME', 0.0  $
                        , 'MODE',  ''   $
                        , 'GRATING', ' ' $
                        , 'AIRMASS', 0.0  $
                        , 'FILEINDX', 0L   $  
                        , 'SKYINDX', -1L $
                        , 'WAVEFRAMES', '')

return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro sinfoni_plan, indir, fileexpr, planfile = planfile, tell_time = tell_time

   COMMON SITE, lat, lng, tzone
   DRADEG = 180.d0/!dpi
   DAY_SEC = 24.0D*60.0D*60.0D
   
   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(fileexpr)) then fileexpr = '*.fits*'
   if (NOT keyword_set(planfile)) then planfile = 'plan.par'
   IF NOT KEYWORD_SET(TELL_TIME) THEN TELL_TIME = 10.0

   spawn, '\ls -d '+indir, dirlist
   if (NOT keyword_set(dirlist)) then begin
       splog, 'No input directories found'
       return
   endif
   ndir = n_elements(dirlist)
   
   ;----------
   ; Loop over each input directory
   observatory, 'eso', obs_struct
   obs_struct.OBSERVATORY =  'par'
   obs_struct.NAME =  'Paranal'
   obs_struct.LONGITUDE = 70.4025d 
   obs_struct.LATITUDE = -24.625d
   obs_struct.ALTITUDE = 2635.43d
   obs_struct.tz = 4.0
   tzone = obs_struct.tz
   lng = 360.d0 - obs_struct.longitude
   lat = obs_struct.latitude
   for idir = 0L, ndir-1L do begin
       splog, 'Working on directory ', dirlist[idir]
       cd, dirlist[idir], current = olddir
       if (idir EQ 0) then origdir = olddir
       filenames = findfile(fileexpr, count = nfile)
       splog, 'Number of FITS files found: ', nfile
       if (nfile GT 0) then begin
          planstr = sinfoni_plan_struct(nfile)
          planstr.filename = filenames
          date_obs = strarr(n_elements(planstr))
          mjd = dblarr(n_elements(planstr))
          ;;expno = lonarr(n_elements(planstr))
          for i = 0L, nfile-1L do begin
              hdr = xheadfits(filenames[i], /silent)
              if (size(hdr, /tname) EQ 'STRING') then begin
                 optic = strcompress(esopar(hdr, 'HIERARCH ESO INS OPTI1 NAME'), /rem)
                 grating = strcompress(esopar(hdr, 'HIERARCH ESO INS GRAT1 NAME'), /rem)
                 exptime = sxpar(hdr, 'EXPTIME')
                 ;;targ = esopar(hdr, 'HIERARCH ESO OBS TARG NAME')
                 targ = strcompress(esopar(hdr, 'HIERARCH ESO OBS TARG NAME'), /rem) $
                        + '-' + strmid(strcompress(esopar(hdr, 'HIERARCH ESO OBS NAME'), /rem), 0, 12)
                 IF targ EQ '-1' THEN targ = $
                    strcompress(esopar(hdr, 'HIERARCH ESO OBS TARG NAME'), /rem)
                 IF targ EQ '-1' THEN targ = strcompress(sxpar(hdr, 'OBJECT'), /rem)
                 catg = strcompress(esopar(hdr, 'HIERARCH ESO DPR CATG'), /rem) ;; science, calib
                 type = strcompress(esopar(hdr, 'HIERARCH ESO DPR TYPE'), /rem) ;; 
                 ;;tech = esopar(hdr, 'HIERARCH ESO DPR TECH')
                 am_sta = esopar(hdr, 'HIERARCH ESO TEL AIRM START')
                 am_end = esopar(hdr, 'HIERARCH ESO TEL AIRM END')
                 airmass = (am_sta + am_end)/2.0d
                 am_str = string((am_sta + am_end)/2.0d, format = '(F5.3)')
                 mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                 planstr[i].SEQ_NUM = long(esopar(hdr, 'HIERARCH ESO TPL EXPNO'))
                 planstr[i].NSEQ = long(esopar(hdr, 'HIERARCH ESO TPL NEXP'))
                 planstr[i].EXPTIME = exptime
                 planstr[i].TARGET = targ
                 planstr[i].GRATING = grating 
                 planstr[i].MODE = optic 
                 planstr[i].AIRMASS  = airmass
                 date_obs[i] = strmid(sxpar(hdr, 'DATE-OBS'), 0, 10)
                 halogen = strcompress(esopar(hdr, 'HIERARCH ESO INS1 LAMP5 ST'), /rem)
                 CASE CATG OF
                    'SCIENCE': BEGIN
                       CASE type OF
                          'OBJECT': planstr[i].FLAVOR = 'science'
                          'SKY':  planstr[i].FLAVOR = 'sky'
                          ELSE: message, 'unknown type'
                       ENDCASE
                    END
                    'CALIB': BEGIN
                       CASE type OF
                          'STD': planstr[i].FLAVOR = 'tell'
                          'SKY,STD':  planstr[i].FLAVOR = 'tell-sky'
                          'FLAT,LAMP': BEGIN
                             CASE halogen OF 
                                'T': planstr[i].FLAVOR = 'iflat-lamp' 
                                'F': planstr[i].FLAVOR = 'iflat-dark'
                                ELSE: message, 'Unknown value for lamp'
                             ENDCASE
                          END
                          'WAVE': planstr[i].FLAVOR = 'arc'
                          'DARK': planstr[i].FLAVOR = 'dark' 
                          ELSE: planstr[i].FLAVOR = 'unknown'
                       ENDCASE
                    END
                    ELSE: planstr[i].FLAVOR = 'unknown'
                 ENDCASE
              ENDIF
           ENDFOR
          isort = sort(mjd)
          planstr = planstr[isort]
          planstr.fileindx = lindgen(nfile)
          mjd = mjd[isort]
          date_obs = date_obs[isort]
          ;;expno = expno[isort]
          ;; Identify different groups and label them 0 to ngroup-1. 
          ;; Grouping is done by mode and grating 
          setup_list = string(planstr.MODE, FORMAT = '(F6.3)') + '-' $
                       + planstr.GRATING ;;+ '-' + date_obs
          isci = WHERE(planstr.FLAVOR EQ 'science' OR planstr.FLAVOR EQ 'tell', nsci1)
          IF nsci1 GT 0 THEN setup_list1 = setup_list[isci] $
          ELSE setup_list1 = setup_list
          setup_uniq = setup_list1[uniq(setup_list1, sort(setup_list1))]
          ngroup = n_elements(setup_uniq)
          FOR j = 0L, ngroup-1L DO BEGIN
             this_group = WHERE(setup_list EQ setup_uniq[j], nthis)
             planstr[this_group].GROUP = j
          ENDFOR
          ;; For each setup, loop over the science frames grouping them into
          ;; offset sequences, and assigning sky and wavelength files
          FOR j = 0L, ngroup-1L DO BEGIN
             group_inds = WHERE(planstr.GROUP EQ j $
                                AND (planstr.FLAVOR EQ 'science' OR  $
                                     planstr.FLAVOR EQ 'sky' OR $
                                     planstr.FLAVOR EQ 'tell' OR $
                                     planstr.FLAVOR EQ 'tell-sky'), nsci)
             IF nsci GT 0 THEN targ_group = planstr[group_inds].TARGET $
             ELSE CONTINUE
             uniq_targ = targ_group[uniq(targ_group)]
             ntarg = n_elements(uniq_targ)
             ;; Loop over every target in this group
             FOR k = 0L, ntarg-1L DO BEGIN
                this_targ = WHERE(planstr.TARGET EQ uniq_targ[k] $
                                  AND planstr.GROUP EQ j $
                                  AND (planstr.FLAVOR EQ 'science' OR  $
                                       planstr.FLAVOR EQ 'sky' OR $
                                       planstr.FLAVOR EQ 'tell' OR $
                                       planstr.FLAVOR EQ 'tell-sky'), nthis)
                nseq = planstr[this_targ[0]].nseq ;; assumes all have same seq
                ;; sort by MJD and then arrange into sequences
                imjd = sort(mjd[this_targ])
                ;;planstr[this_targ[imjd]].SEQ_NUM =lindgen(nthis) MOD nseq
                ;;planstr[this_targ].SEQ_NUM = expno[this_targ] ;; use
                ;;info in header
                iseq = 1
                seq_ctr = 1L
                FOR ll = 0L, nthis-1L DO BEGIN
                   IF planstr[this_targ[imjd[ll]]].SEQ_NUM EQ seq_ctr THEN BEGIN 
                      planstr[this_targ[imjd[ll]]].ISEQ = iseq
                      seq_ctr =  seq_ctr + 1L
                   ENDIF ELSE BEGIN
                      iseq = iseq + 1L ;; increment sequence
                      seq_ctr = 2L     ;; reset counter
                      planstr[this_targ[imjd[ll]]].ISEQ = iseq
                   ENDELSE
                ENDFOR
                ;;planstr[this_targ[imjd]].ISEQ = lindgen(nthis)/nseq
                ;; Check if we have sky frames in this sequence
                isky = WHERE(planstr[this_targ].FLAVOR EQ 'sky' OR $
                             planstr[this_targ].FLAVOR EQ 'tell-sky', nsky)
                ;; Assign sky files for each file of this group/target
                IF nsky GT 0 THEN BEGIN
                   FOR ithis = 0L, nthis-1L DO BEGIN
                      IF planstr[this_targ[ithis]].FLAVOR EQ 'science' OR $
                         planstr[this_targ[ithis]].FLAVOR EQ 'tell' THEN BEGIN
                         diff_mjd = abs(mjd[this_targ[isky]]-mjd[this_targ[ithis]])
                         min_mjd = min(diff_mjd, skyindx)
                         planstr[this_targ[ithis]].SKYINDX = $
                            planstr[this_targ[isky[skyindx]]].FILEINDX
                      ENDIF 
                   ENDFOR
                ENDIF ELSE BEGIN 
                   ;; this uses the sequence design to assign a sky
                   seq_sky = sinfoni_sequence(planstr[this_targ].seq_num, nseq)
                   FOR ithis = 0L, nthis-1L DO BEGIN
                      skyindx = WHERE(planstr[this_targ].seq_num EQ $
                                      seq_sky[ithis] $
                                      AND planstr[this_targ].ISEQ EQ $
                                      planstr[this_targ[ithis]].ISEQ, nsky_ind)
                      IF nsky_ind EQ 0 THEN BEGIN
                         ;; if the whole sequence was not executed, find
                         ;; nearest MJD exluding the image itself
                         diff_mjd = abs(mjd[this_targ]-mjd[this_targ[ithis]])
                         diff_mjd[ithis] = 1d10
                         min_mjd = min(diff_mjd, skyindx)
                      ENDIF ELSE IF nsky_ind GT 1 THEN $
                         message, 'Problem here, nsky_ind should be 1 or zero'
                      planstr[this_targ[ithis]].SKYINDX = $
                         planstr[this_targ[skyindx]].FILEINDX
                   ENDFOR
                ENDELSE
                ;; Assign wavelength files for each science/telluric file of this
                ;; group/target
                this_targ = WHERE(planstr.TARGET EQ uniq_targ[k] AND $
                                  planstr.GROUP EQ j AND $
                                  (planstr.FLAVOR EQ 'science' OR  $
                                   planstr.FLAVOR EQ 'tell'), nthis)
                FOR ithis = 0L, nthis-1L DO BEGIN
                   ;; For 0.25 mode, use object frame for wavelengths
                   IF planstr[this_targ[ithis]].MODE EQ '0.25' THEN BEGIN
                      IF planstr[this_targ[ithis]].FLAVOR EQ 'science' THEN $
                         planstr[this_targ[ithis]].WAVEFRAMES = $
                         sinfoni_frames(planstr[this_targ[ithis]].FILEINDX) $
                      ELSE IF planstr[this_targ[ithis]].FLAVOR EQ 'tell' THEN BEGIN
                         wave_inds = WHERE(planstr.GROUP EQ j $
                                           AND (planstr.FLAVOR EQ 'science'), nwave)
                         IF nwave GT 0 THEN BEGIN 
                            diff_mjd = abs(mjd[this_targ[ithis]]-mjd[wave_inds])
                            min_mjd = min(diff_mjd, kmin)
                            planstr[this_targ[ithis]].WAVEFRAMES = $
                               sinfoni_frames(planstr[wave_inds[kmin]].FILEINDX)
                         ENDIF
                      ENDIF 
                   ENDIF ELSE IF planstr[this_targ[ithis]].MODE EQ '0.025' THEN BEGIN
                      IF planstr[this_targ[ithis]].FLAVOR EQ 'science'  THEN BEGIN 
                         ;; For 0.025 mode, stack all sky frames from the
                         ;; sequence
                         wave_inds = WHERE(planstr.GROUP EQ j AND $
                                           planstr.TARGET EQ uniq_targ[k] AND $
                                           planstr.FLAVOR EQ 'sky' AND $
                                           planstr.ISEQ EQ planstr[this_targ[ithis]].ISEQ, nwave)
                         IF nwave GT 0 THEN planstr[this_targ[ithis]].WAVEFRAMES = $
                            sinfoni_frames(planstr[wave_inds].FILEINDX)
                      ENDIF ELSE IF planstr[this_targ[ithis]].FLAVOR EQ 'tell' THEN BEGIN
                         ;; Grab all possible sky frames with this setup
                         sky_inds = WHERE(planstr.GROUP EQ j AND $
                                           planstr.FLAVOR EQ 'sky', nwave)
                         IF nwave GT 0 THEN BEGIN
                            ;; Find the closest sky frame to the telluric
                            diff_mjd = abs(mjd[this_targ[ithis]]-mjd[sky_inds])
                            min_mjd = min(diff_mjd, kmin)
                            targ_min = planstr[sky_inds[kmin]].TARGET
                            iseq_min = planstr[sky_inds[kmin]].ISEQ
                            wave_inds = WHERE(planstr[sky_inds].TARGET EQ targ_min AND $
                                              planstr[sky_inds].ISEQ EQ iseq_min, nwave)
                            IF nwave GT 0 THEN planstr[this_targ[ithis]].WAVEFRAMES = $
                               sinfoni_frames(planstr[sky_inds[wave_inds]].FILEINDX)
                         ENDIF
                      ENDIF
                   ENDIF ELSE message, 'Current setup not yet supported'
                ENDFOR
             ENDFOR
          ENDFOR
                    
          logfile = repstr(planfile, '.par', '') + '.log'
          plotfile = repstr(planfile, '.par', '') + '.ps'
          
          hdr = ''
          hdr = [hdr, '# Sequences grouped by MODE [0.025, 0.10, 0.25] + GRATING [H,H+K,K]']
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
   
   return
end
;------------------------------------------------------------------------------
