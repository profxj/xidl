;+
; NAME:
;   niri_plan
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



FUNCTION niri_skyframes, inds
skystr = '[' + strcompress(string(inds[0]), /rem)
FOR j = 1L, n_elements(inds)-1L DO BEGIN
    skystr = skystr + ', ' + strcompress(string(inds[j]), /rem)
ENDFOR    
skystr = skystr + ']'

RETURN, skystr
END

FUNCTION niri_plan_struct, nfile
planstr = create_struct(name = 'lexp' $
                        , 'FILENAME', ''  $  
                        , 'GEMINDX', 0L   $  
                        , 'FLAVOR', ''    $
                        , 'TARGET', ''    $
                        , 'EXPTIME', 0.0  $
                        , 'FILTER',  ''   $
                        , 'SLIT',    ''   $
                        , 'GROUP', 0L     $
                        , 'AIRMASS', 0.0  $
                        , 'SEQ_NUM', 0L   $
                        , 'OFFSET', 0.0   $
                        , 'SKYFRAMES', '')

return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro niri_plan, indir, fileexpr, planfile = planfile, tell_time = tell_time

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

   for idir = 0L, ndir-1L do begin
       splog, 'Working on directory ', dirlist[idir]
       cd, dirlist[idir], current = olddir
       if (idir EQ 0) then origdir = olddir
       filenames = findfile(fileexpr, count = nfile)
       splog, 'Number of FITS files found: ', nfile
       
       if (nfile GT 0) then begin
           planstr = niri_plan_struct(nfile)
           planstr.filename = filenames
           mjd_obs = dblarr(n_elements(planstr))
           qoffset = dblarr(n_elements(planstr))
           date_obs = strarr(n_elements(planstr))
           for i = 0L, nfile-1L do begin
               hdr = headfits(filenames[i])
               if (size(hdr, /tname) EQ 'STRING') then begin
                   gcalshut = strtrim(sxpar(hdr, 'GCALSHUT'))
                   filter3 =  strtrim(sxpar(hdr, 'FILTER3'))
                   slit   =  strtrim(sxpar(hdr, 'FPMASK'))
                   obstype = strtrim(sxpar(hdr, 'OBSTYPE'))
                   exptime = float(strtrim(sxpar(hdr, 'EXPTIME')))
                   airmass = float(strtrim(sxpar(hdr, 'AIRMASS')))
                   planstr[i].GEMINDX =  gemindx_extract(filenames[i])
                   planstr[i].EXPTIME = exptime
                   planstr[i].TARGET = strtrim(sxpar(hdr, 'OBJECT'))
                   planstr[i].FILTER = filter3
                   planstr[i].SLIT = slit
                   planstr[i].AIRMASS  = airmass
                   mjd_obs[i]          = double(sxpar(hdr, 'MJD_OBS'))
                   date_obs[i]         = sxpar(hdr, 'DATE-OBS')
                   qoffset[i]          = sxpar(hdr, 'QOFFSET')
                   IF (obstype EQ 'FLAT' AND gcalshut EQ 'CLOSED') OR $
                     (obstype EQ 'DARK') THEN  planstr[i].FLAVOR = 'DARK' $
                   ELSE IF  obstype EQ 'FLAT' AND gcalshut EQ 'OPEN' $
                     THEN  planstr[i].FLAVOR = 'FLAT' $
                   ELSE IF obstype EQ 'ARC'  THEN  planstr[i].FLAVOR = 'ARC'  $
                   ELSE IF (filter3 EQ 'pupil38_G5207' OR slit EQ $
                            'f6-cam_G5208') THEN planstr[i].FLAVOR = 'ACQ' $
                   ELSE IF obstype EQ 'OBJECT' AND  exptime LE TELL_TIME THEN $
                     planstr[i].FLAVOR = 'TELL' $
                   ELSE IF obstype EQ 'OBJECT' AND exptime GE TELL_TIME $
                     THEN planstr[i].FLAVOR = 'SCIENCE' $
                   ELSE message, 'Unrecognized file'
               ENDIF
           ENDFOR
;          Identify different groups and label them 0 to ngroup-1
;          at the moment, grouping is done by date. One could also 
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
;          offset sequences
           FOR j = 0L, ngroup-1L DO BEGIN
               group_inds = WHERE(planstr.GROUP EQ j $
                                  AND (planstr.FLAVOR EQ 'SCIENCE' OR  $
                                       planstr.FLAVOR EQ 'TELL'), nsci)
               IF nsci GT 0 THEN targets = planstr[group_inds].TARGET $
               ELSE CONTINUE
               uniq_targets = targets[uniq(targets, sort(targets))]
               ntarg = n_elements(uniq_targets)
               FOR k = 0L, ntarg-1L DO BEGIN
                   this_targ = WHERE(planstr.TARGET EQ uniq_targets[k] $
                                     AND planstr.GROUP EQ j $ 
                                     AND (planstr.FLAVOR EQ 'SCIENCE' OR  $
                                          planstr.FLAVOR EQ 'TELL'), ntarg)
;                  Sort these images with respect to index
                   this_indx =  planstr[this_targ].GEMINDX
                   indx_sort = sort(this_indx)
                   this_targ = this_targ[indx_sort]
                   this_q    = qoffset[this_targ]
                   planstr[this_targ].OFFSET = this_q-this_q[0]
                   ref_offset = round(planstr[this_targ[0]].OFFSET)
                   seq_beg = WHERE(round(planstr[this_targ].OFFSET) EQ $
                                   ref_offset OR  $
                                   (this_indx-shift(this_indx, 1) GT 1), nseq)
                   seq_end = shift(seq_beg, -1) - 1L 
                   seq_end[n_elements(seq_end)-1L] = n_elements(this_targ)-1L
;                  Loop over the sequences to assign sequence number
                   FOR m = 0L, nseq-1L DO $
                     planstr[this_targ[seq_beg[m]:seq_end[m]]].SEQ_NUM = m
;                  Now loop over the targets strings and assign the 
;                  skyframes for each sequence
                   FOR m = 0L, nseq-1L DO BEGIN
                       this_seq = WHERE(planstr[this_targ].SEQ_NUM EQ m, nthis)
                       indseq = planstr[this_targ[this_seq]].GEMINDX
                       FOR l = 0L, nthis-1L DO BEGIN
                           me_ind = indseq[l]
                           tempo = WHERE(indseq NE me_ind, ntemp)
                           IF ntemp EQ 0 THEN CONTINUE
                           other_ind = indseq[tempo]
                           planstr[this_targ[this_seq[l]]].SKYFRAMES = $
                                   niri_skyframes(other_ind)
                       ENDFOR
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
