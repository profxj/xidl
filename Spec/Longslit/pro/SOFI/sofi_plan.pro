;+
; NAME:
;   sofi_plan
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

FUNCTION sofi_frames, inds
  skystr = '[' + strcompress(string(inds[0]), /rem)
  FOR j = 1L, n_elements(inds)-1L DO BEGIN
     skystr = skystr + ', ' + strcompress(string(inds[j]), /rem)
ENDFOR    
  skystr = skystr + ']'
  
  RETURN, skystr
END


FUNCTION sofi_sequence, seq_num, nseq
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


FUNCTION sofi_plan_struct, nfile
  planstr = create_struct(name = 'lexp'     $
                          , 'FILENAME', ''  $  
                          , 'GROUP', -1L    $
                          , 'SEQ_INDX', 1L      $    
                          , 'STEP_IN_SEQ', 0L   $
                          , 'NSEQ', 0L      $
                          , 'NOD', 0L       $
                          , 'OFFSET', ''    $
                          , 'LOCATION', '-' $
                          , 'FLAVOR', ''    $
                          , 'TARGET', ''    $
                          , 'EXPTIME', 0.0  $
                          , 'MODE',  ''     $
                          , 'SLIT', ''      $
                          , 'AIRMASS', 0.0  $
                          , 'FILEINDX', 0L  $
                          , 'WAVEINDX', 0L  $
                          , 'DARKINDX', 0L)

return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro sofi_plan, indir, fileexpr, planfile = planfile, tell_time = tell_time

   COMMON SITE, lat, lng, tzone
   DRADEG = 180.d0/!dpi
   DAY_SEC = 24.0D*60.0D*60.0D
   
   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(fileexpr)) then fileexpr = '*.fits*'
   if (NOT keyword_set(planfile)) then planfile = 'plan.par'
   IF NOT KEYWORD_SET(TELL_TIME) THEN TELL_TIME = 10.0
   alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

   spawn, '\ls -d '+indir, dirlist
   if (NOT keyword_set(dirlist)) then begin
       splog, 'No input directories found'
       return
   endif
   ndir = n_elements(dirlist)
   
   ;----------
   ; Loop over each input directory
   observatory, 'eso', obs_struct
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
         planstr = sofi_plan_struct(nfile)
         planstr.filename = filenames
         date_obs = strarr(n_elements(planstr))
         mjd = dblarr(n_elements(planstr))
         step_in_seq_ob = lonarr(n_elements(planstr))
         nseq_ob = lonarr(n_elements(planstr))
         ;;expno = lonarr(n_elements(planstr))
          for i = 0L, nfile-1L do begin
             hdr = xheadfits(filenames[i], /silent)
             if (size(hdr, /tname) EQ 'STRING') then begin
                setup = strcompress(esopar(hdr, 'HIERARCH ESO INS MODE'), /rem)
                CASE setup OF
                   'LONG_SLIT_RED': mode = 'H+K'
                   'LONG_SLIT_BLUE': mode = 'J+H'
                   'LONG_SLIT_H': mode = 'H'
                   'LONG_SLIT_K': mode = 'K'
                   'LARGE_FIELD_IMAGING': mode = 'imag'
                   'DARK': mode = 'dark'

                   ELSE: message, 'unrecognized setup'
                ENDCASE
                slit_str1 = strcompress(esopar(hdr, 'HIERARCH ESO INS OPTI1 ID'), /rem)
                IF strmatch(slit_str1, '*large_field*') THEN slit_str = 'open' $
                ELSE slit_str = strmid(slit_str1, 10, 3)
                exptime = sxpar(hdr, 'EXPTIME')
                targ = esopar(hdr, 'HIERARCH ESO OBS TARG NAME')
                name = esopar(hdr, 'HIERARCH ESO OBS NAME')
                IF targ EQ '-1' THEN targ = strcompress(sxpar(hdr, 'OBJECT'), /rem)
                catg = strcompress(esopar(hdr, 'HIERARCH ESO DPR CATG'), /rem) ;; science, calib
                type = strcompress(esopar(hdr, 'HIERARCH ESO DPR TYPE'), /rem) ;;
                dome_str = strcompress(esopar(hdr, 'HIERARCH ESO TEL DOME STATUS'), /rem)
                ;;tech = esopar(hdr, 'HIERARCH ESO DPR TECH')
                am_sta = esopar(hdr, 'HIERARCH ESO TEL AIRM START')
                am_end = esopar(hdr, 'HIERARCH ESO TEL AIRM END')
                airmass = (am_sta + am_end)/2.0d
                am_str = string((am_sta + am_end)/2.0d, format = '(F5.3)')                
                mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                ;; Need this for figuring out if these are twilight flats
                jd = 2400000.5D + mjd[i]
                sunpos, jd, ra1, dec1  ; returns degrees
                zenpos, jd, ra2, dec2  ; returns radians
                ra2 = ra2 * DRADEG
                dec2 = dec2 * DRADEG
                sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                step_in_seq_ob[i] =  long(esopar(hdr, 'HIERARCH ESO TPL EXPNO'))
                planstr[i].STEP_IN_SEQ = step_in_seq_ob[i]
                nseq_ob[i] = long(esopar(hdr, 'HIERARCH ESO TPL NEXP'))
                planstr[i].NSEQ = nseq_ob[i]
                planstr[i].NOD = long(esopar(hdr, 'HIERARCH ESO SEQ NODTHROW'))
                offset = strcompress(esopar(hdr, 'HIERARCH ESO SEQ CUMOFFSETX'), /rem)
                offset = float(offset)
                planstr[i].OFFSET = strcompress(string(offset, format = '(F6.1)'), /rem)
                planstr[i].EXPTIME = exptime
                planstr[i].TARGET = targ
                planstr[i].MODE = mode
                planstr[i].AIRMASS  = airmass
                planstr[i].SLIT  = slit_str
                date_obs[i] = strmid(sxpar(hdr, 'DATE-OBS'), 0, 10)
                halogen = strcompress(esopar(hdr, 'HIERARCH ESO INS1 LAMP5 ST'), /rem)
                CASE CATG OF
                   'SCIENCE': BEGIN
                      ;;(sun_angle GT -12.) THEN $
                      IF (type EQ 'SKY') OR strmatch(name, '*SpecTwi*') THEN $
                         planstr[i].flavor = 'twiflat-sky' $
                      ELSE IF type EQ 'OBJECT' AND exptime LE 20.0 THEN $
                         planstr[i].FLAVOR = 'tell' $
                      ELSE planstr[i].FLAVOR = 'science'
                   END
                   'CALIB': BEGIN
                      CASE type OF
                         'DARK': BEGIN
                            IF (strmatch(name, '*SpecTwi*') OR $
                                dome_str EQ 'FULLY-OPEN') THEN $
                                   planstr[i].FLAVOR = 'twiflat-dark' $
                            ELSE planstr[i].FLAVOR = 'dark'
                         END
                         'FLAT': BEGIN
                            ;; Annoyingly the headers don't say
                            ;; whether the lamp is on or not!!
                            sofi_proc, planstr[i].FILENAME, imag
                            djs_iterstat, imag, median = median
                            IF median GT 1000 THEN planstr[i].FLAVOR = 'dflat-lamp' $
                            ELSE planstr[i].FLAVOR = 'dflat-dark'
                         END
                         'LAMP': planstr[i].FLAVOR = 'arc'
                         'OTHER': BEGIN
                            IF setup EQ 'LARGE_FIELD_IMAGING' AND dome_str EQ 'FULLY-OPEN' THEN $
                               planstr[i].FLAVOR = 'acq' $
                            ELSE planstr[i].FLAVOR = 'other'
                         END  
                         ELSE: message, 'Unrecognized calibration type'
                      ENDCASE
                   END
                   'ACQUISITION': planstr[i].FLAVOR = 'acq'
                   ;; These are typically flats and darks. Not sure
                   ;; why they label them differently
                   'TEST':  planstr[i].FLAVOR = 'test'
                      ;;BEGIN
                      ;IF type EQ 'FLAT' THEN BEGIN
                      ;   sofi_proc, planstr[i].FILENAME, imag
                      ;   djs_iterstat, imag, median = median
                      ;   IF median GT 1000 THEN planstr[i].FLAVOR = 'iflat-lamp' $
                      ;   ELSE planstr[i].FLAVOR = 'iflat-dark'
                      ;ENDIF ELSE planstr[i].FLAVOR = 'test'
                   ;END
                   ELSE: planstr[i].FLAVOR = 'unknown'
                ENDCASE
             ENDIF
          ENDFOR
          isort = sort(mjd)
          planstr = planstr[isort]
          planstr.fileindx = lindgen(nfile)
          mjd = mjd[isort]
          date_obs = date_obs[isort]
          nseq_ob = nseq_ob[isort]
          step_in_seq_ob = step_in_seq_ob[isort]
          ;; Loop over the twilfiat-dark and assign them the setup of
          ;; the nearest twiflat-sky
          idark = where(planstr.FLAVOR EQ 'twiflat-dark', ndark)
          FOR dd = 0L, ndark-1L DO BEGIN
             itwi = WHERE(planstr.FLAVOR EQ 'twiflat-sky' AND $
                          planstr.SLIT EQ planstr[idark[dd]].SLIT AND $
                          planstr.EXPTIME EQ planstr[idark[dd]].EXPTIME, ntwi)
             IF ntwi GT 0 THEN BEGIN 
                diff_mjd = abs(mjd[idark[dd]]-mjd[itwi])
                min_mjd = min(diff_mjd, kmin)
                planstr[idark[dd]].MODE = planstr[itwi[kmin]].MODE
             ENDIF
          ENDFOR
          ;;expno = expno[isort]
          ;; Identify different groups and label them 0 to ngroup-1. 
          ;; Grouping is done by mode and grating 
          setup_list = string(planstr.MODE) 
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
                                     planstr.FLAVOR EQ 'tell'), nsci)
             IF nsci GT 0 THEN targ_group = planstr[group_inds].TARGET $
             ELSE CONTINUE
             uniq_targ = targ_group[uniq(targ_group)]
             ntarg = n_elements(uniq_targ)
             ;; Loop over every target in this group
             FOR k = 0L, ntarg-1L DO BEGIN
                this_targ = WHERE(planstr.TARGET EQ uniq_targ[k] $
                                  AND planstr.GROUP EQ j $
                                  AND (planstr.FLAVOR EQ 'science' OR  $
                                       planstr.FLAVOR EQ 'tell'), nthis)
                ;; For science  frames  8 = 4xA + 4xB sequence
                ;; For telluric frames  2 =  A + B
                this_nod = planstr[this_targ].NOD
                iuniq = uniq(this_nod, sort(this_nod))
                uniq_nod = this_nod[iuniq]
                uniq_nod_mjd = mjd[this_targ[iuniq]]
                uniq_nod = uniq_nod[uniq(uniq_nod_mjd, sort(uniq_nod_mjd))]
                no_of_seq = n_elements(uniq_nod)
                ;; Now loop over the files in each sequence and assign
                ;; nseq and positions
                FOR iseq = 0L, no_of_seq-1L DO BEGIN
                   this_seq = WHERE(planstr[this_targ].NOD EQ uniq_nod[iseq], nseq)
                   mjd_this_seq = mjd[this_targ[this_seq]] ;; MJD sort to set seqeunce numbering
                   this_seq = this_seq[sort(mjd_this_seq)]
                   planstr[this_targ[this_seq]].SEQ_INDX = iseq + 1L ;; assign sequence number
                   planstr[this_targ[this_seq]].NSEQ = nseq          ;; number of steps at this nod
                   planstr[this_targ[this_seq]].STEP_IN_SEQ = lindgen(nseq) + 1L
                   ;; Now identify the offsets in this sequence and
                   ;; assign them a letter A, B (C, D, ... etc.)
                   offsets = planstr[this_targ[this_seq]].OFFSET
                   alpha_indx = 0
                   uniq_offset = offsets[0]
                   FOR ll = 0L, nseq-1L DO BEGIN
                      this_offset = planstr[this_targ[this_seq[ll]]].OFFSET
                      imatch = WHERE(uniq_offset EQ this_offset, nmatch)
                      IF nmatch EQ 0 THEN BEGIN
                         uniq_offset = [uniq_offset, this_offset]
                         alpha_indx = alpha_indx + 1
                         planstr[this_targ[this_seq[ll]]].LOCATION = alphabet[alpha_indx]
                      ENDIF ELSE planstr[this_targ[this_seq[ll]]].LOCATION = alphabet[imatch]
                   ENDFOR
                ENDFOR
                ;;iuniq_offsets = uniq(offsets, sort(offsets))
                ;;uniq_offsets = offsets[iuniq_offsets] 
                ;;uniq_offsets_mjd = mjd[this_targ[this_seq[iuniq_offsets]]]
                ;; MJD sort to set location order
                ;;uniq_offsets = uniq_offsets[uniq(uniq_offsets_mjd, sort(uniq_offsets_mjd))]
                ;;noff = n_elements(uniq_offsets)
                ;;IF planstr[this_targ[0]].FLAVOR EQ 'science' THEN STOP
                ;;FOR ioff = 0L, noff-1L DO BEGIN
                ;;  this_offset = WHERE(planstr[this_targ[this_seq]].OFFSET EQ uniq_offsets[ioff])
                ;;  planstr[this_targ[this_seq[this_offset]]].LOCATION = alphabet[ioff]
                ;;ENDFOR
                ;;ENDFOR
                ;;IF planstr[this_targ[0]].FLAVOR EQ 'science' THEN planstr[this_targ].NSEQ = 8 $
                ;;ELSE planstr[this_targ].NSEQ = 2 ;;nseq_ob[this_targ[0]]
                ;;IF strmatch(planstr[this_targ[0]].TARGET, '*1621-0042_bef*') THEN STOP
                ;;nseq = planstr[this_targ[0]].nseq 
                ;; sort by MJD and then arrange into sequences
                ;;imjd = sort(mjd[this_targ])
                ;planstr[this_targ[imjd]].STEP_IN_SEQ = (lindgen(nthis) MOD nseq) + 1
                ;planstr[this_targ[imjd]].SEQ_INDX = lindgen(nthis)/nseq + 1
                ;;seq_indx = 1
                ;;seq_ctr = 1L
                ;;nod_now = planstr[this_targ[imjd[0]]].NOD
                ;;FOR ll = 0L, nthis-1L DO BEGIN
                ;;   IF planstr[this_targ[imjd[ll]]].NOD EQ nod_now THEN BEGIN
                ;;      planstr[this_targ[imjd[ll]]].STEP_IN_SEQ = seq_ctr
                ;;      planstr[this_targ[imjd[ll]]].SEQ_INDX = seq_indx
                ;;      seq_ctr =  seq_ctr + 1L
                ;;   ENDIF ELSE BEGIN
                ;;      seq_indx = seq_indx + 1L ;; increment sequence
                ;;      planstr[this_targ[imjd[ll]]].SEQ_INDX = seq_indx
                ;;      planstr[this_targ[imjd[ll]]].STEP_IN_SEQ = 1L ;; first throw in new sequence
                ;;      nod_now =  planstr[this_targ[imjd[ll]]].NOD
                ;;      seq_ctr = 2L     ;; reset counter
                ;;   ENDELSE
                ;;ENDFOR
                ;;ainds = where((planstr[this_targ].STEP_IN_SEQ MOD 2) EQ 1, comp = binds)
                ;;planstr[this_targ[ainds]].LOCATION = 'A'
                ;;planstr[this_targ[binds]].LOCATION = 'B'
                ;; 
                ;; Now set the wavelength files for the telluric
                ;; 
                FOR ithis = 0L, nthis-1L DO BEGIN
                   IF planstr[this_targ[ithis]].FLAVOR EQ 'tell' THEN BEGIN
                      ;; telluric uses nearest science frame for
                      ;; wavelength calibration
                      wave_inds = WHERE(planstr.GROUP EQ j $
                                        AND (planstr.FLAVOR EQ 'science'), nwave)
                      IF nwave GT 0 THEN BEGIN 
                         diff_mjd = abs(mjd[this_targ[ithis]]-mjd[wave_inds])
                         min_mjd = min(diff_mjd, kmin)
                         planstr[this_targ[ithis]].WAVEINDX = planstr[wave_inds[kmin]].FILEINDX
                      ENDIF
                   ENDIF ELSE $  ;; science data uses itself for wavelength calibration
                      planstr[this_targ[ithis]].WAVEINDX = planstr[this_targ[ithis]].FILEINDX
                ENDFOR
             ENDFOR             
          ENDFOR
          iflat = WHERE(planstr.FLAVOR EQ 'dflat-lamp', nflat)
          ;; Loop over all the internal flats and assign the dark frames
          FOR kk = 0L, nflat-1L DO BEGIN
             idark = WHERE(planstr.FLAVOR EQ 'dflat-dark' AND $
                           planstr.MODE EQ planstr[iflat[kk]].MODE AND $
                           planstr.SLIT EQ planstr[iflat[kk]].SLIT AND $
                           planstr.EXPTIME EQ planstr[iflat[kk]].EXPTIME, ndark)
             IF ndark GT 0 THEN BEGIN 
                diff_mjd = abs(mjd[iflat[kk]]-mjd[idark])
                min_mjd = min(diff_mjd, kmin)
                planstr[iflat[kk]].DARKINDX = planstr[idark[kmin]].FILEINDX
             ENDIF
          ENDFOR
          iflat = WHERE(planstr.FLAVOR EQ 'twiflat-sky', nflat)
          FOR kk = 0L, nflat-1L DO BEGIN
             idark = WHERE(planstr.FLAVOR EQ 'twiflat-dark' AND $
                           planstr.MODE EQ planstr[iflat[kk]].MODE AND $
                           planstr.SLIT EQ planstr[iflat[kk]].SLIT AND $
                           planstr.EXPTIME EQ planstr[iflat[kk]].EXPTIME, ndark)
             IF ndark GT 0 THEN BEGIN 
                diff_mjd = abs(mjd[iflat[kk]]-mjd[idark])
                min_mjd = min(diff_mjd, kmin)
                planstr[iflat[kk]].DARKINDX = planstr[idark[kmin]].FILEINDX
             ENDIF
          ENDFOR
          logfile = repstr(planfile, '.par', '') + '.log'
          plotfile = repstr(planfile, '.par', '') + '.ps'
          hdr = ''
          hdr = [hdr, '# Sequences grouped by setup [H,K,H+K]']
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
