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

function NIRSPEC_SIGN, A, B
  
  NP   = N_PARAMS()
  npt  = n_elements( A )
  if NP eq 1 then begin
     n    = A
     m    = replicate( 1.d0,npt )
  endif else begin
     n    = B
     m    = A
  endelse
  if npt eq 1 then begin
     if n ge 0 then $
        return, m(0) $
     else           $
        return,-m(0)
  endif else begin
     signs   = A*0
     here_ge = WHERE( n ge 0, nge )
     here_lt = WHERE( n lt 0, nlt )
     if nge gt 0 then signs( here_ge ) =  m( here_ge )
     if nlt gt 0 then signs( here_lt ) = -m( here_lt )
     return, signs
  endelse
end

FUNCTION keckindx_extract, filename
  split1 = strsplit(filename, 's*.fits*', /extract)
  day = strmid(split1[0],3,5) 
  nsplit = n_elements(split1)
  split2 = strmid(split1[nsplit-1L],1,3)
  tmp=day[0] + split2[0]
  RETURN, long(tmp)
END


FUNCTION nirspec_skyframes, inds
skystr = '[' + strcompress(string(inds[0]), /rem)
FOR j = 1L, n_elements(inds)-1L DO BEGIN
    skystr = skystr + ', ' + strcompress(string(inds[j]), /rem)
ENDFOR    
skystr = skystr + ']'

RETURN, skystr
END

FUNCTION nirspec_plan_struct, nfile
planstr = create_struct(name = 'lexp' $
                        , 'FILENAME', ''  $  
                        , 'KECKINDX', 0L   $  
                        , 'FLAVOR', ''    $
                        , 'TARGET', ''    $
                        , 'EXPTIME', 0.0  $
                        , 'FILTER',  ''   $
                        , 'DISPPOS', 0.0  $
                        , 'SLIT',    ''   $
                        , 'GROUP', 0L     $
                        , 'AIRMASS', 0.0  $
                        , 'SEQ_NUM', 0L   $
                        , 'OFFSET', 0.0   $
                        , 'SKYFRAMES', '')

return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro nirspec_plan, indir, fileexpr, planfile = planfile, tell_time = tell_time

   COMMON SITE, lat, lng, tzone
   DRADEG = 180.d0/!dpi
   DAY_SEC = 24.0D*60.0D*60.0D
   
   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(fileexpr)) then fileexpr = '*s*.fits*'
   if (NOT keyword_set(planfile)) then planfile = 'plan.par'
   IF NOT KEYWORD_SET(TELL_TIME) THEN TELL_TIME = 10.0

   spawn, '\ls -d '+indir, dirlist
   if (NOT keyword_set(dirlist)) then begin
       splog, 'No input directories found'
       return
   endif
   ndir = n_elements(dirlist)
   if ndir GT 1 then stop ;; JXP -- I cant get this spawn thing to work..
   dirlist = indir
   
   ;----------
   ; Loop over each input directory
   observatory, 'keck', obs_struct
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
           planstr = nirspec_plan_struct(nfile)
           planstr.filename = filenames
           ra_tel = dblarr(n_elements(planstr))
           dec_tel = dblarr(n_elements(planstr))
           ra_targ = dblarr(n_elements(planstr))
           dec_targ = dblarr(n_elements(planstr))
           date_obs = strarr(n_elements(planstr))
           for i = 0L, nfile-1L do begin
              hdr = headfits(filenames[i],/silent)
              if (size(hdr, /tname) EQ 'STRING') then begin
                 calpos = strtrim(sxpar(hdr, 'CALMPOS'))
                 slit   =  strtrim(sxpar(hdr, 'SLITNAME'))
                 obstype = strtrim(sxpar(hdr, 'COMMENT'))
                 itime = float(strtrim(sxpar(hdr, 'ITIME')))
                 coadds=long(strtrim(sxpar(hdr, 'COADDS')))
                 exptime=itime*float(coadds)
                 flat = long(sxpar(hdr,'FLAT'))
                 neon = long(sxpar(hdr,'NEON'))
                 argon = long(sxpar(hdr,'ARGON'))
                 krypton = long(sxpar(hdr,'KRYPTON'))
                 xenon = long(sxpar(hdr,'XENON'))
                 etalon = long(sxpar(hdr,'ETALON'))
                 am_str = sxpar(hdr,'AIRMASS')
                 ra_tl_str=sxpar(hdr,'RA')
                 dec_tl_str=sxpar(hdr,'DEC')
                 ra_tg_str=sxpar(hdr,'TARGRA')
                 dec_tg_str=sxpar(hdr,'TARGDEC')
                 mjd=string(sxpar(hdr, 'MJD-OBS'))
                 IF strmatch(mjd,'*Error*') THEN sun_angle=-100 $ 
                 ELSE BEGIN
                    jd = 2400000.5D + double(mjd)
                    sunpos, jd, ra1, dec1       ; returns degrees
                    zenpos, jd, ra2, dec2       ; returns radians
                    ra2 = ra2 * DRADEG
                    dec2 = dec2 * DRADEG
                    sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                 ENDELSE
                 
                 IF strmatch(am_str,'*#*') THEN airmass=0.0 $
                 ELSE airmass=float(am_str)
                 planstr[i].KECKINDX = keckindx_extract(filenames[i])
                 planstr[i].EXPTIME = exptime
                 targ =  strtrim(sxpar(hdr, 'TARGNAME'))
                 IF NOT strmatch(targ,'*#*') THEN $
                    planstr[i].TARGET = targ $
                 ELSE planstr[i].TARGET = ''
                 planstr[i].FILTER = strtrim(sxpar(hdr, 'FILNAME'))
                 planstr[i].SLIT = slit
                 planstr[i].AIRMASS  = airmass
                 planstr[i].DISPPOS = float(sxpar(hdr, 'DISPPOS'))
                 date_obs[i]         = sxpar(hdr, 'DATE-OBS')
                 IF NOT strmatch(ra_tl_str,'*#*') THEN ra_tel[i]=ra_tl_str
                 IF NOT strmatch(dec_tl_str,'*#*') THEN dec_tel[i]=dec_tl_str
                 IF NOT strmatch(ra_tg_str,'*#*') THEN ra_targ[i]=ra_tg_str
                 IF NOT strmatch(dec_tg_str,'*#*') THEN dec_targ[i]=dec_tg_str
                 IF strmatch(planstr[i].FILTER,'*BLANK*') THEN BEGIN
                    IF exptime LE 5 THEN planstr[i].FLAVOR = 'dark' $
                    ELSE planstr[i].FLAVOR = 'longdark'
                 ENDIF ELSE IF (calpos EQ 1) THEN BEGIN
                    IF FLAT THEN planstr[i].FLAVOR = 'iflat' $
                    ELSE IF NEON THEN planstr[i].FLAVOR = 'arc-ne' $
                    ELSE IF ARGON THEN planstr[i].FLAVOR = 'arc-ar' $
                    ELSE IF KRYPTON THEN planstr[i].FLAVOR = 'arc-kr' $
                    ELSE IF XENON THEN planstr[i].FLAVOR = 'arc-xe'  $
                    ELSE IF ETALON THEN planstr[i].FLAVOR = 'arc-et' $
                    ELSE IF exptime LE 5 THEN planstr[i].FLAVOR = 'dark' $
                    ELSE planstr[i].FLAVOR = 'unknown'
                 ENDIF ELSE BEGIN
                    IF (sun_angle GT -12.) then $
                        planstr[i].flavor = 'twiflat' $
                    ELSE IF exptime LE TELL_TIME THEN $
                       planstr[i].FLAVOR = 'tell' $
                    ELSE planstr[i].FLAVOR = 'science'
                 ENDELSE
              ENDIF
           ENDFOR
;          Identify different groups and label them 0 to ngroup-1. 
;          Grouping is done by cross-disperser angle. 
           setup_list = string(planstr.DISPPOS, FORMAT = '(F6.3)') + '-' $
                        + planstr.SLIT + '-' + date_obs
           setup_uniq = setup_list[uniq(setup_list,sort(setup_list))]
           ;;date_obs, sort(date_obs))]
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
                                      planstr.FLAVOR EQ 'tell') $
                                 AND NOT (ra_targ EQ 0.0 AND dec_targ EQ 0.0) $
                                 , nsci)
              IF nsci GT 0 THEN ra_targ_group = ra_targ[group_inds] $
              ELSE CONTINUE
              uniq_ra_targ = ra_targ_group[uniq(ra_targ_group $
                                                , sort(ra_targ_group))]
              ntarg = n_elements(uniq_ra_targ)
              FOR k = 0L, ntarg-1L DO BEGIN
                 this_targ = WHERE(ra_targ EQ uniq_ra_targ[k] $
                                   AND planstr.GROUP EQ j $ 
                                   AND (planstr.FLAVOR EQ 'science' OR  $
                                        planstr.FLAVOR EQ 'tell'), nthis)
;                Sort these images with respect to index
                 this_indx =  planstr[this_targ].KECKINDX
                 indx_sort = sort(this_indx)
                 this_targ = this_targ[indx_sort] 
;                 IF (this_indx[sort(this_indx)])(0) EQ 7099 THEN STOP
                 this_ra=ra_tel[this_targ]
                 this_dec=dec_tel[this_targ]
                 aoffset=3600.0*djs_diff_angle(this_ra,this_dec $
                                               ,this_ra[0],this_dec[0])
                 dec_rad=dec_targ[this_targ[0]]*!dpi/180.0
                 cos_dec = cos(dec_rad)
                 ra_off=3600.0*cos_dec*(this_ra - this_ra[0])
                 dec_off=3600.0*(this_dec-this_dec[0])
                 IF total(abs(ra_off))/double(nthis) GT 0.2 THEN $
                    sgn=nirspec_sign(ra_off) $
                 ELSE sgn = nirspec_sign(dec_off)
                 planstr[this_targ].OFFSET = sgn*aoffset
                 ref_offset = round(planstr[this_targ[0]].OFFSET)
                 seq_beg = WHERE(round(planstr[this_targ].OFFSET) EQ $
                                 ref_offset OR  $
                                 (this_indx-shift(this_indx, 1) GT 1), nseq)
;                 IF nthis/nseq GT 5 THEN BEGIN
;                    nseq = ceil(float(nthis)/5)
;                    nper=nthis/nseq
;                    seq_beg=[seq_beg[0],seq_beg[0] + (lindgen(nseq))[1:*]*nper]
;                 ENDIF
                 seq_end = shift(seq_beg, -1) - 1L 
                 seq_end[n_elements(seq_end)-1L] = n_elements(this_targ)-1L
;                Loop over the sequences to assign sequence number
                 FOR m = 0L, nseq-1L DO $
                    planstr[this_targ[seq_beg[m]:seq_end[m]]].SEQ_NUM = m
;                  Now loop over the targets strings and assign the 
;                  skyframes for each sequence
                   FOR m = 0L, nseq-1L DO BEGIN
                       this_seq = WHERE(planstr[this_targ].SEQ_NUM EQ m, nthis)
                       indseq = planstr[this_targ[this_seq]].KECKINDX
                       FOR l = 0L, nthis-1L DO BEGIN
                          me_ind = indseq[l]
                          tempo = WHERE(indseq NE me_ind, ntemp)
                          IF ntemp EQ 0 THEN CONTINUE
                          other_ind = indseq[tempo]
                          diff=min(abs(other_ind-me_ind),kmin)
                          ioth=WHERE(abs(other_ind-me_ind) EQ diff,noth)
                          IF noth GT 1 THEN $
                             ioth=WHERE(abs(other_ind-me_ind) EQ diff AND $
                                        other_ind-me_ind LT 0)
                          planstr[this_targ[this_seq[l]]].SKYFRAMES = $
                             nirspec_skyframes(other_ind[ioth])
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
