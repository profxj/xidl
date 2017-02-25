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
;   xheadfits()
;   idlutils_version()
;   splog
;   sxpar()
;   yanny_write
;
; INTERNAL SUPPORT ROUTINES:
;   niri_plan_struct()
;
; REVISION HISTORY:
;   27-Jun-2008  Written by Joseph Hennawi, Berkeley
;-
;------------------------------------------------------------------------------


FUNCTION mage_plan_struct, nfile
planstr = create_struct(name = 'lexp' $
                        , 'FILENAME', ''  $  
                        , 'FLAVOR', ''    $
                        , 'TARGET', ''    $
                        , 'EXPTIME', 0.0  $
                        , 'SLIT',    ''   $
                        , 'GROUP', 0L     $
                        , 'INSTRUMENT', ''$
                        , 'AIRMASS', 0.0)
return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro mage_plan, indir, fileexpr, planfile = planfile, tell_time = tell_time

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
         planstr = mage_plan_struct(nfile)
         planstr.filename = filenames
         mjd_obs = dblarr(n_elements(planstr))
         date_obs = strarr(n_elements(planstr))
         for i = 0L, nfile-1L do begin
               hdr = xheadfits(filenames[i])
               if (size(hdr, /tname) EQ 'STRING') then begin
                  exptime = sxpar(hdr, 'EXPTIME')
                  date_obs = sxpar(hdr, 'DATE-OBS')
                  ut_time  = sxpar(hdr, 'UT-TIME')
                  ;; Decide if this is twilight: sun_angle is the
                  ;; angle of the sun *above* the horizon, so must
                  ;; be less than -12 if darker than 12-degree twi.
                  observatory, 'lco', obs_struct
                  tzone = obs_struct.tz
                  lng = 360.d0 - obs_struct.longitude
                  lat = obs_struct.latitude
                  jd = x_setjdate(date_obs, ut_time)
                  sunpos, jd, ra1, dec1      ; returns degrees
                  zenpos, jd, ra2, dec2      ; returns radians
                  ra2 = ra2 * DRADEG
                  dec2 = dec2 * DRADEG
                  sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                  exp_type = strtrim(sxpar(hdr, 'EXPTYPE'))
                  if strmatch(exp_type, '*Flat*') THEN $
                     planstr[i].flavor = 'domeflat' $
                  else if strmatch(exp_type, '*ThAr-Lamp*') THEN $
                     planstr[i].flavor = 'arc' $
                  else if strmatch(exp_type, '*Xe-Flash*') THEN $
                     planstr[i].flavor = 'flashflat' $
                  else if (sun_angle GT -12.) then $
                     planstr[i].flavor = 'twiflat' $
                  else if exptime LT 120.0 THEN $
                     planstr[i].flavor = 'std' $
                  else if strmatch(exp_type, '*Object*') THEN $
                     planstr[i].flavor = 'science' $
                  ELSE message, 'Unidentified file type'
                  planstr[i].SLIT = strtrim(sxpar(hdr, 'SLITNAME'))
                  planstr[i].exptime = exptime
                  planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                  instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                  planstr[i].instrument = instrument
                  airmass = float(strtrim(sxpar(hdr, 'AIRMASS')))
                  planstr[i].AIRMASS  = airmass
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
