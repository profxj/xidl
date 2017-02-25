;+
; NAME:
;   luci_plan
;
; PURPOSE:
;   Create plan file(s) for running the LUCI pipeline
;
; CALLING SEQUENCE:
;   luci_plan, [ fileexpr, indir, planfile= ]
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

FUNCTION luci_sequence, seq_num, nseq
  CASE nseq OF
     ;; Hack to deal with cases where NSEQ not set in headers
     0: sky_ind = replicate(-1, n_elements(seq_num))
     1: BEGIN 
        ;; calibration image
        sky_ind = [0] 
     END
     2: BEGIN 
        ;; 2-pt Dither pattern, each uses the other
        sky_ind = [1, 0] 
     END
     4: BEGIN 
        ;; 4-pt Dither pattern
        ;;  2     0     3     1
        ;;  none uses 0 image as sky. 2 and 3 use each other
        ;;         0, 1, 2, 3 
        sky_ind = [2, 3, 4, 3] 
     END
     6: BEGIN 
        ;; 6-pt Dither pattern
        ;; 
        ;;  none uses 0 image as sky. 2 and 3 use each other
        ;;         0, 1, 2, 3, 4, 5, 6
        sky_ind = [1, 2, 3, 4, 5, 6, 5]
     END  
     8: BEGIN 
        ;;  Rockabily sequences are labeled as -1
        ;;sky_ind = [1, 2, 3, 4, 5, 6, 7, 6]
        sky_ind = replicate(-1, 8)
     END  
     ELSE: message, 'Unsupported sequence number. Add it here'
  ENDCASE
  RETURN, sky_ind[seq_num-1]
END


FUNCTION luci_indx_extract, files

  nfiles = n_elements(files)
  indx = lonarr(nfiles)
  FOR ii = 0L, nfiles-1L DO BEGIN 
     split1 = strsplit(files[ii], 'luci*.fits*', /extract)
     IF n_elements(split1) GT 1 THEN tmp = long(split1[1]) $
     ELSE tmp = 0
     indx[ii] = long(tmp)
  ENDFOR
  RETURN, indx
END


FUNCTION luci_plan_struct, nfile
  planstr = create_struct(name = 'lexp' $
                          , 'FILENAME', ''  $  
                          , 'GROUP', -1L     $
                          , 'SEQ_NUM', 0L   $
                          , 'ISEQ', 0L   $    
                          , 'NSEQ', 0L   $  
                          , 'FLAVOR', ''    $
                          , 'TARGET', ''    $
                          , 'EXPTIME', 0.0  $
                          , 'GRATING',  ''   $
                          , 'WCEN', 0.0  $
                          , 'SLIT',    ''   $
                          , 'AIRMASS', 0.0  $
                          , 'FILEINDX', 0L   $  
                          , 'SKYINDX', 0L)
  
return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro luci_plan, indir, fileexpr, planfile = planfile, tell_time = tell_time
  
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
   observatory, 'mgio', obs_struct
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
          planstr = luci_plan_struct(nfile)
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
                 slit_str = strcompress(sxpar(hdr, 'MASKID'), /rem)
                 CASE slit_str OF
                    '990065': slit = 'slit_0.25'
                    '990078': slit = 'slit_0.50'
                    '990029': slit = 'slit_0.75'
                    '990034': slit = 'slit_1.00'
                    'LS0.5_300mue': slit = 'slit_0.50'
                    'LS_0.50arcsec': slit = 'slit_0.50'
                    'blind': slit = 'blind'
                    ELSE: message, 'ERROR: Unknown mask'
                 ENDCASE
                 exptime = sxpar(hdr, 'ITIME')
                 targ = sxpar(hdr, 'OBJECT')
                 ;;tech = esopar(hdr, 'HIERARCH ESO DPR TECH')
                 grat = strcompress(sxpar(hdr, 'GRATNAME'), /rem)
                 wcen = strcompress(string(sxpar(hdr, 'GRATWLEN'), format = '(F5.3)'), /rem)
                 ;; Compute the airmass of the observation
                 ra_deg = double(sxpar(hdr, 'CRVAL1'))
                 dec_deg = double(sxpar(hdr, 'CRVAL2'))
                 lst_sec = double(sxpar(hdr, 'LST'))
                 lst_hr = lst_sec/(60.0d*60.0d)
                 hour_angle = ra_deg/15.0d - lst_hr
                 airm = airmass(lat, dec_deg, hour_angle)
                 am_str = string(airm, format = '(F5.3)')
                 mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                 jd = 2400000.5D + mjd[i]
                 sunpos, jd, ra1, dec1 ; returns degrees
                 zenpos, jd, ra2, dec2 ; returns radians
                 ra2 = ra2 * DRADEG
                 dec2 = dec2 * DRADEG
                 sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                 ;; ????? LUCI2 commisioning doesn't have these
                 ;; ????? header cards for NSEQ? Check to see if LUCI2 final headers
                 ;; ???? do? 
                 planstr[i].NSEQ = long(sxpar(hdr, 'NELEM'))
                 planstr[i].SEQ_NUM = long(sxpar(hdr, 'ELEM_NBR'))
                 planstr[i].EXPTIME = exptime
                 planstr[i].TARGET = targ
                 ;;planstr[i].FILTER = strtrim(sxpar(hdr, 'FILNAME'))
                 planstr[i].SLIT = slit
                 planstr[i].AIRMASS  = airm
                 planstr[i].GRATING = grat 
                 planstr[i].WCEN = wcen
                 date_obs[i] = strmid(sxpar(hdr, 'DATE-OBS'), 0, 10)
                 ne_lamp = strmatch(sxpar(hdr, 'STATLMP1'), '*ON*')
                 ar_lamp = strmatch(sxpar(hdr, 'STATLMP2'), '*ON*')
                 xe_lamp = strmatch(sxpar(hdr, 'STATLMP3'), '*ON*')
                 halo1_lamp = strmatch(sxpar(hdr, 'STATLMP4'), '*ON*')
                 halo2_lamp = strmatch(sxpar(hdr, 'STATLMP5'), '*ON*')
                 halo3_lamp = strmatch(sxpar(hdr, 'STATLMP6'), '*ON*')
                 filter1 = strcompress(sxpar(hdr, 'FILTER1'), /rem)
                 filter2 = strcompress(sxpar(hdr, 'FILTER2'), /rem)
                 telmode = strcompress(sxpar(hdr, 'TELMODE'), /rem)
                 readmode = strcompress(sxpar(hdr, 'READMODE'), /rem)
                 nexpo = long(sxpar(hdr, 'NEXPO'))
                 IF ne_lamp OR ar_lamp OR xe_lamp THEN $
                    planstr[i].flavor = 'arc' $
                 ELSE IF halo1_lamp OR halo2_lamp OR halo3_lamp THEN $
                    planstr[i].flavor = 'iflat-lamp' $
                 ELSE IF strmatch(filter1, '*blind*') OR strmatch(filter2, '*blind*') THEN $
                    planstr[i].flavor = 'dark' $
                 ELSE IF strmatch(grat, '*mirror*') THEN BEGIN
                    IF strmatch(sxpar(hdr, 'MOSPOS'), '*Turnout*') THEN $
                       planstr[i].FLAVOR = 'acq-imag' $
                    ELSE IF strmatch(sxpar(hdr, 'MOSPOS'), '*FPU*') THEN $
                       planstr[i].FLAVOR = 'acq-slit' $
                    ELSE planstr[i].FLAVOR = 'acq-unknown' 
                 ENDIF ELSE IF strmatch(TELMODE, '*TRACK*') THEN BEGIN
                    ;; Could be a dark or a twilight flat?
                    IF (sun_angle GT -8. AND sun_angle LT 8.0) THEN $
                       planstr[i].flavor = 'twiflat' $
                    ELSE IF nexpo GT 1 THEN planstr[i].FLAVOR = 'iflat-dark' $
                    ELSE planstr[i].FLAVOR = 'tell'
                 ENDIF ELSE BEGIN
                    IF strmatch(readmode, '*o2.double.corr.read*') THEN $
                       planstr[i].FLAVOR = 'tell' $
                    ELSE planstr[i].FLAVOR = 'science'
                 ENDELSE
              ENDIF
           ENDFOR           
           isort = sort(mjd)
           planstr = planstr[isort]
           planstr.fileindx = luci_indx_extract(planstr.FILENAME)
           mjd = mjd[isort]
           date_obs = date_obs[isort]
;          Identify different groups and label them 0 to ngroup-1. 
;          Grouping is done by cross-disperser angle. 
           setup_list = planstr.GRATING + '-' + $
                        string(planstr.WCEN, FORMAT = '(F6.3)') + '-' $
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
                 ;;planstr[this_targ[imjd]].SEQ_NUM = lindgen(nthis) MOD nseq
                 planstr[this_targ[imjd]].ISEQ = lindgen(nthis)/nseq 
                 ;; this uses the sequence design to assign a sky
                 seq_sky = luci_sequence(planstr[this_targ].seq_num, nseq)
                 ;; Is this a rockabilly sequence?
                 dum = where(seq_sky EQ -1, nrock)
                 IF nrock EQ n_elements(seq_sky) THEN $
                    planstr[this_targ].SKYINDX = -1L $
                 ELSE BEGIN 
                    FOR ithis = 0L, nthis-1L DO BEGIN
                       ;;IF strmatch(planstr[this_targ[0]].TARGET, '*HIP*') THEN STOP
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
                 ENDELSE
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
