;+
; NAME:
;   isaac_plan
;
; PURPOSE:
;   Create plan file(s) for running the GNIRS pipeline
;
; CALLING SEQUENCE:
;   niri_plan, [ fileexpr, indir, planfile= ]
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

FUNCTION isaac_sequence, seq_num, nseq
  CASE nseq OF
     2: BEGIN 
        ;; 2-pt Dither pattern, each uses the other
        sky_ind = [1, 0] 
     END
     4: BEGIN 
        ;; 4-pt Dither pattern
        ;;  2     0     3     1
        ;;  none uses 0 image as sky. 2 and 3 use each other
        sky_ind = [1, 2, 3, 2] 
     END
     6: BEGIN 
        ;; 6-pt Dither pattern
        ;; 
        ;;  none uses 0 image as sky. 2 and 3 use each other
        sky_ind = [1, 2, 3, 4, 5, 6, 5]
     END  
     8: BEGIN 
        ;;  none uses 0 image as sky. 2 and 3 use each other
        sky_ind = [1, 2, 3, 4, 5, 6, 7, 6] 
     END  
     ELSE: message, 'Unsupported sequence number. Add it here'
  ENDCASE
  RETURN, sky_ind[seq_num]
END


FUNCTION nirspec_skyframes, inds
skystr = '[' + strcompress(string(inds[0]), /rem)
FOR j = 1L, n_elements(inds)-1L DO BEGIN
    skystr = skystr + ', ' + strcompress(string(inds[j]), /rem)
ENDFOR    
skystr = skystr + ']'

RETURN, skystr
END

FUNCTION isaac_plan_struct, nfile
planstr = create_struct(name = 'lexp' $
                        , 'FILENAME', ''  $  
                        , 'GROUP', -1L     $
                        , 'SEQ_NUM', 0L   $
                        , 'ISEQ', 0L   $    
                        , 'NSEQ', 0L   $  
                        , 'FLAVOR', ''    $
                        , 'TARGET', ''    $
                        , 'EXPTIME', 0.0  $
                        , 'MODE',  ''   $
                        , 'WCEN', 0.0  $
                        , 'SLIT',    ''   $
                        , 'AIRMASS', 0.0  $
                        , 'FILEINDX', 0L   $  
                        , 'SKYINDX', 0L)

return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro isaac_plan, indir, fileexpr, planfile = planfile, tell_time = tell_time

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
           planstr = isaac_plan_struct(nfile)
           planstr.filename = filenames
           ra_tel = dblarr(n_elements(planstr))
           dec_tel = dblarr(n_elements(planstr))
           ra_targ = dblarr(n_elements(planstr))
           dec_targ = dblarr(n_elements(planstr))
           date_obs = strarr(n_elements(planstr))
           mjd = dblarr(n_elements(planstr))
           for i = 0L, nfile-1L do begin
              hdr = xheadfits(filenames[i], /silent)
              if (size(hdr, /tname) EQ 'STRING') then begin
                 slit_str = $
                    strcompress(esopar(hdr, 'HIERARCH ESO INS OPTI1 ID'), /rem)
                 CASE slit_str OF 
                    'slit_2': slit = 'slit_2.0'
                    'slit_1': slit = 'slit_1.0'
                    'slit_0.6_tilted': slit = 'slit_0.6'
                    ELSE: slit = strcompress(slit_str, /rem)
                 ENDCASE
                 exptime = sxpar(hdr, 'EXPTIME')
                 targ = esopar(hdr, 'HIERARCH ESO OBS TARG NAME')
                 IF targ EQ '-1' THEN targ = sxpar(hdr, 'OBJECT')
                 catg = esopar(hdr, 'HIERARCH ESO DPR CATG') ;; science, calib
                 type = esopar(hdr, 'HIERARCH ESO DPR TYPE') ;; 
                 ;;tech = esopar(hdr, 'HIERARCH ESO DPR TECH')
                 grat = esopar(hdr, 'HIERARCH ESO INS GRAT NAME')
                 wcen = strcompress( $
                        string(esopar(hdr, 'HIERARCH ESO INS GRAT WLEN') $
                               , format = '(F5.3)'), /rem)
                 order = strcompress(long(esopar(hdr, 'HIERARCH ESO INS GRAT ORDER')), /rem)
                 calshut = esopar(hdr, 'HIERARCH ESO INS CALSHUT ST')
                 CASE order OF 
                    '2': band = 'K'
                    '3': band = 'H'
                    '4': band = 'J'
                    ELSE: message, 'unrecognized order'
                 ENDCASE
                 am_sta = esopar(hdr, 'HIERARCH ESO TEL AIRM START')
                 am_end = esopar(hdr, 'HIERARCH ESO TEL AIRM END')
                 airmass = (am_sta + am_end)/2.0d
                 am_str = string((am_sta + am_end)/2.0d, format = '(F5.3)')
                 mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                 jd = 2400000.5D + mjd[i]
                 sunpos, jd, ra1, dec1    ; returns degrees
                 zenpos, jd, ra2, dec2    ; returns radians
                 ra2 = ra2 * DRADEG
                 dec2 = dec2 * DRADEG
                 sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                 ;;planstr[i].SEQEXPNO = long(esopar(hdr, 'HIERARCH ESO TPL EXPNO'))
                 planstr[i].NSEQ = long(esopar(hdr, 'HIERARCH ESO TPL NEXP'))
                 ;;planstr[i].KECKINDX = keckindx_extract(filenames[i])
                 planstr[i].EXPTIME = exptime
                 planstr[i].TARGET = targ
                 ;;planstr[i].FILTER = strtrim(sxpar(hdr, 'FILNAME'))
                 planstr[i].SLIT = slit
                 planstr[i].AIRMASS  = airmass
                 planstr[i].MODE = strcompress(grat, /rem)  + '-' + band
                 planstr[i].WCEN = wcen
                 date_obs[i] = strmid(sxpar(hdr, 'DATE-OBS'), 0, 10)
                 IF (strmatch(type, '*OBJECT*') OR $ 
                     strmatch(catg, '*SCIENCE*')) AND $
                    (sun_angle GT -12.) THEN $
                       planstr[i].flavor = 'twiflat' $
                 ELSE IF (strmatch(type, '*OBJECT*') OR $ 
                          strmatch(catg, '*SCIENCE*')) THEN $
                             planstr[i].FLAVOR = 'science' $
                 ELSE IF strmatch(type, '*STD*') THEN $
                    planstr[i].FLAVOR = 'tell' $
                 ELSE IF strmatch(type, '*FLAT*') THEN BEGIN
                    IF calshut EQ 'T' THEN $ 
                       planstr[i].FLAVOR = 'iflat,lamp' $
                    ELSE IF calshut EQ 'F' THEN $
                       planstr[i].FLAVOR = 'iflat,dark' $
                    ELSE planstr[i].FLAVOR = 'iflat'
                 ENDIF ELSE IF strmatch(type, '*WAVE*') THEN $
                    planstr[i].FLAVOR = 'arc'  $
                 ELSE IF strmatch(type, '*DARK*') THEN $
                    planstr[i].FLAVOR = 'dark' $
                 ELSE IF strmatch(type, '*STARTRACE*') THEN $
                    planstr[i].FLAVOR = 'eso-startrace' $
                 ELSE IF strmatch(type, '*SLIT*') THEN $
                    planstr[i].FLAVOR = 'acq-slit' $
                 ELSE IF strmatch(type, '*SKY*') THEN $
                    planstr[i].FLAVOR = 'acq-imag' $
                 ELSE planstr[i].FLAVOR = 'unknown'
              ENDIF
           ENDFOR
           isort = sort(mjd)
           planstr = planstr[isort]
           planstr.fileindx = lindgen(nfile)
           mjd = mjd[isort]
           date_obs = date_obs[isort]
;          Identify different groups and label them 0 to ngroup-1. 
;          Grouping is done by cross-disperser angle. 
           setup_list = string(planstr.WCEN, FORMAT = '(F6.3)') + '-' $
                        + planstr.SLIT + '-' + date_obs
           isci = WHERE(planstr.FLAVOR EQ 'science', nsci1)
           IF nsci1 GT 0 THEN setup_list1 = setup_list[isci] $
           ELSE setup_list1 = setup_list
           setup_uniq = setup_list1[uniq(setup_list1, sort(setup_list1))]
           ngroup = n_elements(setup_uniq)
           FOR j = 0L, ngroup-1L DO BEGIN
              this_group = WHERE(setup_list EQ setup_uniq[j], nthis)
              planstr[this_group].GROUP = j
           ENDFOR
;          For each setup, loop over the science frames grouping them into
;          offset sequences
           FOR j = 0L, ngroup-1L DO BEGIN
              group_inds = WHERE(planstr.GROUP EQ j $
                                 AND (planstr.FLAVOR EQ 'science' OR  $
                                      planstr.FLAVOR EQ 'tell'), nsci)
              IF nsci GT 0 THEN targ_group = planstr[group_inds].TARGET $
              ELSE CONTINUE
              uniq_targ = targ_group[uniq(targ_group)]
              ntarg = n_elements(uniq_targ)
              FOR k = 0L, ntarg-1L DO BEGIN
                 this_targ = WHERE(planstr.TARGET EQ uniq_targ[k] $
                                   AND planstr.GROUP EQ j $ 
                                   AND (planstr.FLAVOR EQ 'science' OR  $
                                        planstr.FLAVOR EQ 'tell'), nthis)
                 nseq = planstr[this_targ[0]].nseq ;; assumes all have same seq
                 ;; sort by MJD and then arrange into sequence
                 imjd = sort(mjd[this_targ])
                 planstr[this_targ[imjd]].SEQ_NUM = lindgen(nthis) MOD nseq
                 planstr[this_targ[imjd]].ISEQ = lindgen(nthis)/nseq 
                 ;; this uses the sequence design to assign a sky
                 seq_sky = isaac_sequence(planstr[this_targ].seq_num, nseq)
                 FOR ithis = 0L, nthis-1L DO BEGIN
                    skyindx = WHERE(planstr[this_targ].seq_num EQ $
                                    seq_sky[ithis] $
                                    AND planstr[this_targ].ISEQ EQ $
                                    planstr[this_targ[ithis]].ISEQ, nsky)
                    IF nsky EQ 0 THEN BEGIN
                       ;; if the whole sequence was not executed, find
                       ;; nearest MJD exluding the image itself
                       diff_mjd = abs(mjd[this_targ]-mjd[this_targ[ithis]])
                       diff_mjd[ithis] = 1d10
                       min_mjd = min(diff_mjd, skyindx)
                    ENDIF ELSE IF nsky GT 1 THEN $
                       message, 'Problem here, nsky should be 1 or zero'
                    planstr[this_targ[ithis]].SKYINDX = $
                       planstr[this_targ[skyindx]].FILEINDX
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
   
   return
end
;------------------------------------------------------------------------------
