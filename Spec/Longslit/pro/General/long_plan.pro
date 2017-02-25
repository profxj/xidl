;+
; NAME:
;   long_plan
;
; PURPOSE:
;   Create plan file(s) for running the low-redux pipeline.  This code
;   parses headers, does image stats, etc.
;
; CALLING SEQUENCE:
;   long_plan, [ fileexpr, indir, planfile= ]
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
;     bias
;     domeflat
;     iflat (internal flat)
;     twiflat
;     arc
;     science
;
; EXAMPLES:
; long_plan,'*.fits','/b/martell/data_arx/09072005/Raw/',planfile='plan-master.par'
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
;   long_plan_struct()
;
;  To Do: for twilight images taken with a grism in but no mask,
;  assume they are slitless, and treat them as pixel flat field images 
;  not twiflats. 
;
;
;
; REVISION HISTORY:
;   13-Mar-2005  Written by David Schlegel, LBL.
;-
;------------------------------------------------------------------------------
function long_plan_struct, nfile
   planstr = create_struct(name='lexp', $
    'FILENAME'    , '', $
    'FLAVOR'      , '', $
    'TARGET'      , '', $
    'EXPTIME'     , 0., $
    'INSTRUMENT'  , '', $
    'GRATING'     , '', $
    'WAVE'        , '', $                       
    'MASKNAME'    , ''  )
   return, replicate(planstr, nfile)
end
;------------------------------------------------------------------------------
pro long_plan, fileexpr, indir, planfile=planfile

   COMMON SITE, lat, lng, tzone
   DRADEG = 180.d0/!dpi
   
   if (NOT keyword_set(indir)) then indir = '.'
   if (NOT keyword_set(fileexpr)) then fileexpr = '*.fits*'
   if (NOT keyword_set(planfile)) then planfile = 'plan.par'

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

   for idir = 0L, ndir-1L do begin
       splog, 'Working on directory ', dirlist[idir]
       cd, dirlist[idir], current = olddir
       if (idir EQ 0) then origdir = olddir
       filenames = findfile(fileexpr, count = nfile)
       splog, 'Number of FITS files found: ', nfile
       
       if (nfile GT 0) then begin
           planstr = long_plan_struct(nfile)
           planstr.filename = filenames
           mjd = dblarr(n_elements(planstr))
           for i = 0L, nfile-1L do begin
              ;print, 'Reading ', filenames[i]
               hdr = xheadfits(filenames[i])
               if (size(hdr, /tname) EQ 'STRING') then begin
                  if (strmatch(sxpar(hdr, 'INSTRUME'), 'LRIS*')) then begin
                      ;;---------------------------
                      ;; This is the LRIS on Keck I
                      ;;---------------------------
                       exptime = sxpar(hdr, 'ELAPTIME')
                       trapdoor = strtrim(sxpar(hdr, 'TRAPDOOR'))
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       ;; Figure out if telescope is pointing
                       ;; throught the dome slit (to within 4
                       ;; degrees), or at the
                       ;; inside of the dome.  Note that domeflats are
                       ;; generally taken with a telescope
                       ;; position 90 degrees from the slit,
                       ;; but we wont force this case to be true.
                       slitposn = sxpar(hdr, 'AZ')
                       domeposn = sxpar(hdr, 'DOMEPOSN')
                       slitposn = (slitposn + 360.) mod 360.
                       domeposn = (domeposn + 360.) mod 360.
                       if (abs(domeposn - slitposn) lt 4.) then $
                          throughslit = 1 else throughslit = 0
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       observatory, 'keck', obs_struct
                       tzone = obs_struct.tz
                       lng = 360.d0 - obs_struct.longitude
                       lat = obs_struct.latitude
                       mjd[i] = sxpar(hdr, 'MJD-OBS')
                       jd = 2400000.5D + mjd[i]
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                       lamps = sxpar(hdr, 'LAMPS')
                       MERCURY = strmatch(sxpar(hdr, 'MERCURY'), '*on*')
                       NEON = strmatch(sxpar(hdr, 'NEON'), '*on*')
                       ARGON = strmatch(sxpar(hdr, 'ARGON'), '*on*')
                       CADMIUM = strmatch(sxpar(hdr, 'CADMIUM'), '*on*')
                       ZINC = strmatch(sxpar(hdr, 'ZINC'), '*on*')
                       HALOGEN = strmatch(sxpar(hdr, 'HALOGEN'), '*on*')
                       KRYPTON = strmatch(sxpar(hdr, 'KRYPTON'), '*on*')
                       XENON = strmatch(sxpar(hdr, 'XENON'), '*on*')
                       FEARGON = strmatch(sxpar(hdr, 'FEARGON'), '*on*')
                       DEUTERI = strmatch(sxpar(hdr, 'DEUTERI'), '*on*')
                       arclamps = MERCURY OR NEON OR ARGON OR CADMIUM OR ZINC OR $
                                  KRYPTON OR XENON OR FEARGON
                       ilamps = DEUTERI OR HALOGEN
                       if (KEYWORD_SET(lamps) AND lamps EQ '0,0,0,0,0,1') OR $
                          ilamps then planstr[i].flavor = 'iflat' $
                       else if (KEYWORD_SET(lamps) AND $
                                strmid(lamps, 0, 9) NE '0,0,0,0,0') OR $
                          arclamps then planstr[i].flavor = 'arc' $
                       else if (trapdoor EQ 'closed') then $
                          planstr[i].flavor = 'bias' $
                       else if (throughslit eq 0) then $
                          planstr[i].flavor = 'domeflat' $
                       else if (sun_angle GT -12.) then $
                         planstr[i].flavor = 'twiflat' $
                       else if exptime LT 120.0 THEN $
                         planstr[i].flavor = 'std' $
                       else planstr[i].flavor = 'science'
                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'TARGNAME'))
                       planstr[i].maskname = strtrim(sxpar(hdr, 'SLITNAME'))
                  ; Note that GRISNAME will only be defined for blue CCDs,
                  ; and GRANAME will only be defined for red CCDs. 
                  ; MSWAVE is only defined for red CCDs
                       instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                       planstr[i].instrument = instrument
                       planstr[i].grating = $
                         + (instrument EQ 'LRISBLUE' ? sxpar(hdr, 'GRISNAME') : '') $
                         + (instrument EQ 'LRIS' ? sxpar(hdr, 'GRANAME') : '')
                       planstr[i].wave = $
                         + (instrument EQ 'LRISBLUE' ? '' : '') $
                         + (instrument EQ 'LRIS' ? $
                            strcompress(sxpar(hdr, 'MSWAVE'), /rem): '')
                       
                    endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*CA-3.5*') $
                    then begin
                       ;;---------------------------
                       ;; This is the CAHA 3.5m
                       ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       observatory, 'ca', obs_struct
                       tzone = obs_struct.tz
                       lng = 360.d0 - obs_struct.longitude
                       lat = obs_struct.latitude
                       jd = sxpar(hdr, 'JUL-DATE')
                       mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                       caha_hh = caha_hdr(hdr)
                       imagetyp = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                       IF imagetyp EQ 'dark' OR imagetyp EQ 'bias' $
                       THEN planstr[i].flavor = 'bias' $
                       ELSE IF imagetyp EQ 'simulation' OR imagetyp EQ 'arc'THEN $
                          planstr[i].flavor = 'arc' $
                       ELSE IF imagetyp EQ 'flat' THEN $
                          planstr[i].flavor = 'domeflat' $
                       ELSE IF imagetyp EQ 'science' THEN BEGIN 
                          IF sun_angle GT -12. THEN $
                             planstr[i].flavor = 'twiflat' $
                          ELSE IF exptime LT 350.0 THEN $
                             planstr[i].flavor = 'std' $
                          ELSE planstr[i].flavor = 'science'
                       ENDIF ELSE message,  'Unrecognized obstype'
                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                       planstr[i].maskname = strtrim(sxpar(caha_hh, 'SLIT_WID'))
                       instr1 = strcompress(sxpar(hdr, 'INSTRUME'), /rem)
                       path = strcompress(sxpar(caha_hh, 'PATH'), /rem)
                       instrument = instr1 + '-' + path
                       planstr[i].instrument = instrument
                       planstr[i].grating = $
                          + (path EQ 'BLUE' ? sxpar(caha_hh, 'GRAT1_NA') : '') $
                          + (path EQ 'RED'  ? sxpar(caha_hh, 'GRAT2_NA') : '')
                       planstr[i].wave = $
                          + (path EQ 'BLUE' ? sxpar(caha_hh, 'GRAT1_WL') : '') $
                          + (path EQ 'RED'  ? sxpar(caha_hh, 'GRAT2_WL') : '')
                    endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*CA-2.2*') $
                    then begin
                       ;;---------------------------
                       ;; This is the CAHA 2.2m
                       ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       observatory, 'ca', obs_struct
                       tzone = obs_struct.tz
                       lng = 360.d0 - obs_struct.longitude
                       lat = obs_struct.latitude
                       jd = sxpar(hdr, 'JUL-DATE')
                       mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                       ;; This may be unnecessary
                       caha_hh = hdr
                       imagetyp = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                       IF imagetyp EQ 'dark' OR imagetyp EQ 'bias' $
                       THEN planstr[i].flavor = 'bias' $
                       ELSE IF imagetyp EQ 'simulation' OR imagetyp EQ 'arc'THEN $
                          planstr[i].flavor = 'arc' $
                       ELSE IF imagetyp EQ 'flat' THEN $
                          planstr[i].flavor = 'domeflat' $
                       ELSE IF imagetyp EQ 'science' THEN BEGIN 
                          IF sun_angle GT -12. THEN $
                             planstr[i].flavor = 'twiflat' $
                          ELSE IF exptime LT 350.0 THEN $
                             planstr[i].flavor = 'std' $
                          ELSE planstr[i].flavor = 'science'
                       ENDIF ELSE message,  'Unrecognized obstype'
                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                       planstr[i].maskname = strtrim(sxpar(caha_hh, 'INSAPDY')) ;; Microns not arcsec
                       instr1 = strcompress(sxpar(hdr, 'INSTRUME'), /rem)
                       path = strcompress(sxpar(caha_hh, 'PATH'), /rem)
                       instrument = instr1 + '-' + path
                       planstr[i].instrument = instrument
                       planstr[i].grating = strtrim(sxpar(caha_hh,'INSGRNAM'),2)
                       planstr[i].wave = float(sxpar(caha_hh,'INSGRWL0'))
                    endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*ESO-NTT*') $
                    then begin
                       ;;---------------------------
                       ;; This is EFOSC on the NTT 3.6m
                       ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       observatory, 'eso', obs_struct
                       tzone = obs_struct.tz
                       lng = 360.d0 - obs_struct.longitude
                       lat = obs_struct.latitude
                       mjd[i] = sxpar(hdr, 'MJD-OBS') 
                       jd = 2400000.5D + mjd[i]
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                       ind_catg = $
                          WHERE(stregex(hdr, 'HIERARCH ESO DPR CATG', /bool))
                       catg = strcompress(repstr(strmid(hdr[ind_catg], 30, 14), "'", ''), /rem)
                       ind_type = $
                          WHERE(stregex(hdr, 'HIERARCH ESO DPR TYPE', /bool))
                       type = strcompress(repstr(strmid(hdr[ind_type], 30, 14), "'", ''), /rem)
                       IF catg EQ 'CALIB' THEN BEGIN
                          IF type EQ 'DARK' OR type EQ 'BIAS' $
                          THEN planstr[i].flavor = 'bias' $
                          ELSE IF type EQ 'WAVE' THEN $
                             planstr[i].flavor = 'arc' $
                          ELSE IF type EQ 'FLAT' THEN $
                             planstr[i].flavor = 'domeflat' 
                       ENDIF ELSE IF catg EQ 'SCIENCE' THEN BEGIN
                          IF sun_angle GT -12. THEN $
                             planstr[i].flavor = 'twiflat' $
                          ELSE IF exptime LT 350.0 OR type EQ 'STD' THEN $
                             planstr[i].flavor = 'std' $
                          ELSE planstr[i].flavor = 'science'
                       ENDIF ELSE message,  'Unrecognized obstype'
                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                       ind_slit = $
                          WHERE(stregex(hdr, 'HIERARCH ESO INS SLIT1 NAME' $
                                        , /bool))
                       slit = strcompress(repstr(strmid(hdr[ind_slit], 30, 14), "'", ''), /rem)
                       slit = repstr(slit, '#', '')
                       planstr[i].maskname = slit
                       planstr[i].instrument = $
                          strcompress(sxpar(hdr, 'INSTRUME'), /rem)
                       ind_gris = $
                          WHERE(stregex(hdr, 'HIERARCH ESO INS GRIS1 NAME' $
                                        , /bool))
                       gris = strcompress(repstr(strmid(hdr[ind_gris], 30, 14), "'", ''), /rem)
                       gris = repstr(gris, '#', '')
                       planstr[i].grating = gris
                       planstr[i].wave = '0.0'
           endif else if  (strmatch(sxpar(hdr, 'INSTRUME'), 'GMOS-N*')) OR $
                     (strmatch(sxpar(hdr, 'INSTRUME'), 'GMOS-S*')) $
                     then begin
                    ;---------------------------
                    ; This is GMOS on Gemini-N/S
                    ;---------------------------
                       exptime = strcompress(sxpar(hdr, 'EXPOSURE'), /rem)
                       obstype = strcompress(sxpar(hdr, 'OBSTYPE'), /rem)
                       object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                       instrument = strcompress(sxpar(hdr, 'INSTRUME'), /rem)
                       mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                       if obstype EQ 'FLAT' then planstr[i].flavor = 'domeflat' $
                       else if obstype EQ 'ARC' then planstr[i].flavor = 'arc' $
                       else if obstype EQ 'BIAS'then planstr[i].flavor = 'bias' $
                       else if (obstype EQ 'OBJECT' AND object EQ 'Twilight') $
                         then planstr[i].flavor = 'twiflat' $
                       else if obstype EQ 'MASK'then planstr[i].flavor $
                          = 'mask' $
                       else if obstype EQ 'DARK'then planstr[i].flavor $
                          = 'dark' $
                       else if (obstype EQ 'OBJECT' AND object NE 'Twilight') $
                         then BEGIN
                           if exptime LT 120.0 THEN planstr[i].flavor = 'std' $
                           ELSE  planstr[i].flavor = 'science'
                       ENDIF ELSE message,  'Unrecognized obstype'
                       planstr[i].exptime = exptime
                       planstr[i].target = object
                       planstr[i].maskname = strtrim(sxpar(hdr, 'MASKNAME'))
                       instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                  planstr[i].instrument = instrument
                  ;; XGMOS ?
                  if sxpar(hdr,'XGMOS') then $
                    planstr[i].instrument = planstr[i].instrument+'X'
                  planstr[i].grating =  strcompress(sxpar(hdr, 'GRATING') $
                                                    , /rem)
                  planstr[i].wave = strcompress(sxpar(hdr, 'GRWLEN'), /rem)
              endif else if (strmatch(sxpar(hdr, 'INSTRUME'), 'DEIMOS*')) then begin
                  ;---------------------------
                  ; This is the DEIMOS on Keck I
                  ;---------------------------
                       exptime = sxpar(hdr, 'ELAPTIME')
                       trapdoor = strtrim(sxpar(hdr, 'TRAPDOOR'))
                       mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       observatory, 'keck', obs_struct
                       tzone = obs_struct.tz
                       lng = 360.d0 - obs_struct.longitude
                       lat = obs_struct.latitude
                       mjd[i] =  sxpar(hdr, 'MJD-OBS')
                       jd = 2400000.5D + mjd[i]
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)

                       lamps = sxpar(hdr, 'LAMPS')
                       if strmatch(lamps,'Qz*') then $
                         planstr[i].flavor = 'domeflat' $
                       else if (strmid(lamps, 0, 3) NE 'Off') then $
                         planstr[i].flavor = 'arc' $
                       else if strmatch(trapdoor,'closed*') then $
                         planstr[i].flavor = 'bias' $
                       else if (sun_angle GT -12.) then $
                         planstr[i].flavor = 'twiflat' $
                       else begin
                           if exptime LT 120.0 THEN $
                             planstr[i].flavor = 'std' $
                           else planstr[i].flavor = 'science'
                       endelse 
                       planstr[i].exptime = exptime
                       if not strmatch(sxpar(hdr,'MOSMODE'),'Spectral') then $
                         planstr[i].flavor = 'unknown'
                       
                       planstr[i].target = strtrim(sxpar(hdr, 'TARGNAME'))
                       planstr[i].maskname = strtrim(sxpar(hdr, 'SLMSKNAM'))

                       instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                       planstr[i].instrument = instrument
                       planstr[i].grating = strtrim(sxpar(hdr, 'GRATENAM'),2)
                       gpos = strtrim(sxpar(hdr, 'GRATEPOS'),2)

                       planstr[i].wave = sxpar(hdr,'G'+strtrim(gpos,2)+'TLTWAV')
              endif else if strmatch(sxpar(hdr, 'TELESCOP'), 'mmt*') THEN BEGIN
                    ;---------------------------
                    ; This is Blue/Red Channel on MMT 
                    ;---------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  obstype = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  ccd    =  strcompress(sxpar(hdr, 'DETECTOR'), /rem)
                  mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                  if ccd EQ 'ccd35' or ccd EQ 'mmtbluechan' then begin
                     instrument = 'MMT Blue Channel' 
                     dispers = sxpar(hdr, 'DISPERSE')
                     IF KEYWORD_SET(DISPERS) THEN $
                        planstr[i].grating = strcompress(dispers, /rem) $
                     ELSE planstr[i].grating =  '831 2nd order'
                     ;planstr[i].wave = ''
                     planstr[i].wave = strtrim(sxpar(hdr, 'CENWAVE'),2)
                     planstr[i].maskname = strtrim(sxpar(hdr, 'APERTURE'))
                  ENDIF else if (ccd EQ 'ccd34') OR strmatch(ccd,'mmtredchan') then BEGIN
                     instrument = 'MMT Red Channel' 
                     ;; These next lines might not work with the old
                     ;; chip
                     dispers = sxpar(hdr, 'DISPERSE')
                     prs = strsplit(dispers,'-',/extract)
                     planstr[i].grating =  prs[0]
                     ;planstr[i].grating = dispers
                     planstr[i].wave = prs[1]
                   
                  Endif ELSE message, 'Unrecognized MMT instrument'
                  if obstype EQ 'flat' then planstr[i].flavor = 'domeflat' $
                  else if obstype EQ 'comp' then planstr[i].flavor = 'arc' $
                  else if obstype EQ 'zero'then planstr[i].flavor = 'bias' $
                  else if obstype EQ 'object' then begin
                      if exptime LT 120.0 THEN planstr[i].flavor = 'std' $
                      else planstr[i].flavor = 'science' 
                  endif else message,  'Unrecognized obstype'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  ;planstr[i].maskname = object
                  planstr[i].instrument = instrument
                  ;stop
               endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*kp4m*') $
               THEN BEGIN
                  ;;---------------------------
                  ;; This is RC Spectrograph on the KPNO 4m
                  ;;--------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  obstype = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)                  
                  ccd    =  strcompress(sxpar(hdr, 'DETECTOR'), /rem)
                  mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                  instrument = 'RC Spectrograph'
                  dispers = sxpar(hdr, 'DISPERSE')
                  planstr[i].grating = strcompress(dispers, /rem) 
                  planstr[i].wave = strcompress(sxpar(hdr, 'TILTPOS'), /rem)
                  if obstype EQ 'flat' then planstr[i].flavor = 'domeflat' $
                  else if obstype EQ 'comp' then planstr[i].flavor = 'arc' $
                  else if obstype EQ 'zero'then planstr[i].flavor = 'bias' $
                  else if obstype EQ 'object' then begin
                     if exptime LT 200.0 THEN planstr[i].flavor = 'std' $
                     else planstr[i].flavor = 'science' 
                  endif else if obstype EQ 'focus' THEN $
                     planstr[i].flavor = 'focus' $
                  else message,  'Unrecognized obstype'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  planstr[i].maskname =  $
                     strcompress(sxpar(hdr, 'APERTURE'), /rem)
                  planstr[i].instrument = instrument
               endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*DuPont*') $
               then begin
                  ;;---------------------------
                  ;; This is the DuPont 100" with the B&C
                  ;;---------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  obstype = strcompress(sxpar(hdr, 'EXPTYPE'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                  instrument = 'B&C'
                  ;maskname = strcompress(sxpar(hdr, 'SLITMASK'), /rem)
                  planstr[i].grating = 'GRATING'
                  planstr[i].wave = ''
                  if strlen(sxpar(hdr,'LAMPS')) GT 2 then planstr[i].flavor = 'arc'  $
                  else begin
                     if strcompress(sxpar(hdr,'EXPTYPE'),/rem) EQ 'Flat' then planstr[i].flavor = 'bothflat' $
                     else begin
                        IF exptime GT 10.0 THEN planstr[i].flavor = 'science' $
                        ELSE  planstr[i].flavor = '????'
                     endelse
                  endelse
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  ;planstr[i].maskname = maskname
                  planstr[i].instrument = instrument
               endif else if strmatch(sxpar(hdr, 'INSTRUME'), 'IMACS*') $
               THEN BEGIN
                    ;---------------------------
                    ; This is IMACS on Magellan 
                    ;---------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  ;obstype = strcompress(sxpar(hdr, 'EXPTYPE'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  instrument = 'IMACS'
                  mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                  maskname = strcompress(sxpar(hdr, 'SLITMASK'), /rem)
                  planstr[i].grating = strcompress(sxpar(hdr, 'DISPERSR'), /rem)
                  planstr[i].wave = strtrim(sxpar(hdr, 'G-ANGLE'), 2)
                  IF exptime GT 100.0 THEN planstr[i].flavor = 'science' $
                  ELSE  planstr[i].flavor = 'calib'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  planstr[i].maskname = maskname
                  planstr[i].instrument = instrument
              endif else if  strmatch(sxpar(hdr, 'DETECTOR'), 'mars*') $
                then begin
                 mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  obstype = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  instrument = 'KPNO Mars'
                  if obstype EQ 'flat' then planstr[i].flavor = 'domeflat' $
                  else if obstype EQ 'comp' then planstr[i].flavor = 'arc' $
                  else if obstype EQ 'zero'then planstr[i].flavor = 'bias' $
                  else if obstype EQ 'object' then $
                  planstr[i].flavor = 'science' $
                  else message,  'Unrecognized obstype'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  planstr[i].maskname = object
                  planstr[i].instrument = instrument
                  planstr[i].grating =  '8050'
                  planstr[i].wave = ''
;              endif else if (strmatch(strtrim(sxpar(hdr, 'INSTRUME'),2), 'Kast') OR strmatch (strtrim(sxpar(hdr, 'INSTRUME'),2), 'Kast blue arm') OR  strmatch(strtrim(sxpar(hdr, 'INSTRUME'),2), 'KAST blue') OR strmatch(strtrim(sxpar(hdr, 'INSTRUME'),2), 'KAST'))  then begin

              endif else if (stregex(sxpar(hdr,'INSTRUME'),'.*kast.*',$
                                    /boolean,/fold_case) eq 1) or $
                (stregex(sxpar(hdr,'VERSION'),'kast*', /bool, /fold_case) EQ 1) then begin
                  ;;---------------------------
                  ;; This is Kast on Shane 3m (Mt. Hamilton)
                  ;;---------------------------
                  mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                  data = xmrdfits(filenames[i],/fscale)
                  exptime = float(sxpar(hdr, 'EXPOSURE')) > $
                            float(sxpar(hdr, 'EXPTIME'))
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  planstr[i].target = object
                  planstr[i].exptime = exptime
                  ;; Blue or Red?
                  if strtrim(sxpar(hdr,'SPSIDE'),2) EQ 'blue' OR $
                    strmid(sxpar(hdr,'VERSION'),0,5) EQ 'kastb' then begin
                      planstr[i].instrument = 'KAST-B' 
                      planstr[i].grating = strtrim(sxpar(hdr,'GRISM_N'),2)
                      planstr[i].wave = sxpar(hdr,'GRISM_P')
                  endif else begin
                      planstr[i].instrument = 'KAST-R' 
                      planstr[i].grating = strtrim(sxpar(hdr,'GRATNG_N'),2)
                      planstr[i].wave = sxpar(hdr,'GRTILT_P')
                  endelse
                  ;; Type
                  if strmatch(planstr[i].grating,'open') then begin
                      print, 'Imaging ', filenames[i], '.  Skipping..'
                      planstr[i].flavor = 'IMG'
                      continue
                  endif
                  if strtrim(sxpar(hdr,'SHUTTER'),2) NE 'OPEN' AND $
                    strtrim(sxpar(hdr,'OBSTYPE'),2) NE 'OBJECT' then begin
                      if exptime LE 2 then planstr[i].flavor = 'bias' $
                      else planstr[i].flavor = 'dark'
                  endif else begin
                      if exptime LT 120. then begin
                          if median(data) GT 500. then planstr[i].flavor = 'domeflat' $
                          else planstr[i].flavor = 'arc'  
                      endif else planstr[i].flavor = 'science'  
                  endelse
              endif else if strmatch(strtrim(sxpar(hdr, 'INSTRUME'), 2) $
                                     , 'DIS') then begin
                  ;---------------------------
                  ; This is DIS at APO
                  ;---------------------------
                 mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJNAME'), /rem)
                  IF strmatch(object, '*TestSlew*') THEN $
                    planstr[i].target = 'calib' $
                  ELSE planstr[i].target = object
                  planstr[i].exptime = exptime
                  ;; Blue or Red?
                  if strmatch(sxpar(hdr, 'DETECTOR'), '*blue*') then $
                    planstr[i].instrument = 'DIS-B' $
                  else if strmatch(sxpar(hdr, 'DETECTOR'), '*red*') then $
                    planstr[i].instrument = 'DIS-R' $
                  else planstr[i].instrument = ' '
                  planstr[i].grating = strtrim(sxpar(hdr, 'GRATING'), 2)
                  planstr[i].wave = sxpar(hdr, 'DISPWC')
                  planstr[i].maskname = strtrim(sxpar(hdr, 'SLITMASK'))
                  if strmatch(sxpar(hdr, 'IMAGETYP'), '*zero*')  $
                    then planstr[i].flavor = 'bias' $
                  else if  strmatch(sxpar(hdr, 'IMAGETYP'), '*flat*') $
                    then planstr[i].flavor = 'iflat' $
                  else if strmatch(sxpar(hdr, 'IMAGETYP'), '*object*') AND $
                    exptime GT 120.0 THEN planstr[i].flavor = 'science' $
                  else if strmatch(sxpar(hdr, 'IMAGETYP'), '*object*') AND $
                    exptime LE 120.0 THEN planstr[i].flavor = 'arc' $
                  else planstr[i].flavor = 'unknown'  
               endif else if strmatch(strtrim(sxpar(hdr, 'INSTRUME'), 2) $
                                      , 'ISIS*') then begin
                  ;;-------------------------------
                  ;; This is ISIS at WHT, La Palma
                  ;; Added by MF on Nov 2015 
                  ;;-------------------------------
                  mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  planstr[i].target = object
                  planstr[i].exptime = exptime
                  ;; Blue or Red?
                  if strmatch(sxpar(hdr, 'ISIARM'), '*Blue*') then $
                    planstr[i].instrument = 'ISIS-B' $
                  else if strmatch(sxpar(hdr, 'ISIARM'), '*Red*') then $
                    planstr[i].instrument = 'ISIS-R' $
                  else planstr[i].instrument = ' '
                  planstr[i].grating = strtrim(sxpar(hdr, 'ISIGRAT'), 2)
                  planstr[i].wave = sxpar(hdr, 'CENWAVE')
                  planstr[i].maskname = strtrim(sxpar(hdr, 'ISISLITW'))
                  if strmatch(sxpar(hdr, 'IMAGETYP'), '*zero*')  $
                  then planstr[i].flavor = 'bias' $
                  else if  strmatch(sxpar(hdr, 'IMAGETYP'), '*flat*') $
                  then planstr[i].flavor = 'iflat' $
                  else if strmatch(sxpar(hdr, 'IMAGETYP'), '*arc*') $
                  then planstr[i].flavor = 'arc' $
                  else if strmatch(sxpar(hdr, 'IMAGETYP'), '*object*') $
                  then planstr[i].flavor = 'science' $
                  else planstr[i].flavor = 'unknown'  
               endif else if strmatch(strtrim(sxpar(hdr, 'TELID'), 2), '200') $
                  THEN BEGIN
                    ;---------------------------
                    ; This is DBSP on Palomar 200"
                    ;---------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  planstr[i].exptime = exptime
                  mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                  ;; Blue or Red?
                  IF strtrim(sxpar(hdr, 'FPA'), 2) EQ 'DBSP_BLUE' THEN $
                     planstr[i].instrument = 'P200-B' $
                  ELSE IF strtrim(sxpar(hdr, 'FPA'), 2) EQ 'DBSP_RED' THEN $
                     planstr[i].instrument = 'P200-R' $
                  ELSE planstr[i].instrument = ' '
                  planstr[i].grating = strtrim(sxpar(hdr, 'GRATING'), 2)
                  planstr[i].wave =  strtrim(sxpar(hdr, 'ANGLE'), 2)
                  planstr[i].maskname = 'slit' + strtrim(sxpar(hdr, 'APERTURE'))
                  lamps  = strtrim(sxpar(hdr, 'LAMPS'), 2)
                  lamp_vec = long(strmid(lamps, lindgen(7), 1))
                  turret = strtrim(sxpar(hdr, 'TURRET'), 2)
                  imgtype=strtrim(sxpar(hdr,'IMGTYPE'))
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  IF strmatch(imgtype, 'flat') THEN BEGIN
                     planstr[i].flavor = 'iflat' 
                     planstr[i].target = object
                  ENDIF ELSE IF strmatch(imgtype, 'cal') THEN BEGIN
                     planstr[i].flavor = 'arc'
                     lamp_str = ''
                     lamp_names = ['D', 'FeAr', 'Hg', 'Ar', 'Ne', 'He'$
                                   , 'InCand']
                     FOR j = 0L, 6 DO BEGIN
                        IF lamp_vec[j] GT 0 THEN BEGIN
                           IF NOT KEYWORD_SET(lamp_str) THEN $
                              lamp_str = lamp_names[j] $
                           ELSE lamp_str = lamp_str + '-' + lamp_names[j]
                        ENDIF
                     ENDFOR
                     planstr[i].target = lamp_str
                  ENDIF ELSE IF exptime EQ 0 THEN BEGIN
                     planstr[i].flavor = 'bias' 
                     planstr[i].target = 'none'
                  ENDIF ELSE IF strmatch(imgtype, 'object') THEN BEGIN
                     if exptime LT 120.0 THEN planstr[i].flavor = 'std' $
                     else planstr[i].flavor = 'science' 
                     planstr[i].target = object
                  ENDIF ELSE planstr[i].flavor = 'unknown'  
               endif else if (strmatch(sxpar(hdr, 'INSTRUME'), '*FORS2*')) then begin
                  ;;---------------------------
                  ;; This is FORS2 on VLT
                  ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       object = strcompress(sxpar(hdr, 'OBJECT'), /rem)
                       targ = esopar(hdr, 'HIERARCH ESO OBS TARG NAME')
                       IF targ EQ '-1' THEN targ = sxpar(hdr, 'OBJECT')
                       catg = esopar(hdr, 'HIERARCH ESO DPR CATG') ;; science, calib
                       type = esopar(hdr, 'HIERARCH ESO DPR TYPE') ;; 
                 ;;tech = esopar(hdr, 'HIERARCH ESO DPR TECH')
                       gris = esopar(hdr, 'HIERARCH ESO INS GRIS1 NAME')
                       ind_chip = WHERE(stregex(hdr, 'EXTNAME', /bool))
                       chipno = strcompress(strmid(hdr[ind_chip],11,5))
                       IF gris EQ '-1' THEN gris = 'none'
                       wcen = strcompress( $
                        string(esopar(hdr, 'HIERARCH ESO INS GRIS1 WLEN') $
                               , format = '(F5.1)'), /rem)
                       IF wcen EQ '-1' THEN wcen = '--' 
                       mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                       ;; determine if this is a twilight exposure
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
                       jd = 2400000.5D + mjd[i]
                       sunpos, jd, ra1, dec1 ; returns degrees
                       zenpos, jd, ra2, dec2 ; returns radians
                       ra2 = ra2 * DRADEG
                       dec2 = dec2 * DRADEG
                       sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                       IF (strmatch(type, '*OBJECT*') OR $ 
                           strmatch(catg, '*SCIENCE*')) AND $
                          (sun_angle GT -12.) THEN $
                             planstr[i].flavor = 'twiflat' $
                       ELSE IF (strmatch(type, '*OBJECT*') OR $ 
                                strmatch(catg, '*SCIENCE*')) THEN $
                                   planstr[i].FLAVOR = 'science' $
                       ELSE IF strmatch(type, '*STD*') THEN $
                          planstr[i].FLAVOR = 'std' $
                       ELSE IF strmatch(type, '*FLAT*') THEN planstr[i].flavor = 'domeflat' $
                       ELSE IF strmatch(type, '*WAVE*') THEN planstr[i].FLAVOR = 'arc'  $
                       ELSE IF strmatch(type, '*DARK*') THEN planstr[i].FLAVOR = 'dark' $
                       ELSE IF strmatch(type, '*BIAS*') THEN planstr[i].FLAVOR = 'bias' $
                       ELSE IF strmatch(type, '*SKY*') THEN planstr[i].FLAVOR = 'acq-imag' $
                       ELSE IF strmatch(type, '*SLIT*') THEN planstr[i].FLAVOR = 'acq-slit' $
                       ELSE planstr[i].FLAVOR = 'unknown'
                       planstr[i].exptime = exptime
                       planstr[i].target = targ 
                       ;; maskname
                       mode = esopar(hdr, 'HIERARCH ESO INS MODE')
                       IF strmatch(mode, '*MXU*') THEN BEGIN
                          mask = esopar(hdr, 'HIERARCH ESO INS MASK NAME') 
                          IF mask EQ '-1' THEN maskname = 'mask-acq' $
                          ELSE maskname = mask
                       ENDIF ELSE IF strmatch(mode, '*MOS*') THEN BEGIN
                          wid = esopar(hdr, 'HIERARCH ESO INS MOS1 WID  ')
                          slit = strcompress(string(wid, format = '(F4.1)'), /rem)
                          pos = esopar(hdr, 'HIERARCH ESO INS MOS1 POS')
                          pos_str = strcompress(string(pos, format = '(I)'), /rem)
                          maskname = 'MOS-slit=' + slit + '-pos=' + pos_str 
                       ENDIF ELSE IF  strmatch(mode, '*MOS*') THEN BEGIN
                          maskname = 'IMG'
                       ENDIF ELSE IF strmatch(mode, '*LSS*') THEN $
                          maskname = 'LSS'
                       planstr[i].maskname = maskname
                       instrument = strtrim(sxpar(hdr, 'INSTRUME'))+chipno
                       planstr[i].instrument = instrument
                       planstr[i].grating = gris
                       planstr[i].wave = wcen
                       
               endif else if (strmatch(sxpar(hdr, 'INSTRUME'), '*MODS1*')) then begin
                      ;;---------------------------
                      ;; This is MODS on LBT
                      ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       imtype = strtrim(sxpar(hdr, 'IMAGETYP'))
                       mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                       if imtype eq 'BIAS' then planstr[i].flavor = 'bias' $
                       else if imtype eq 'COMP' then planstr[i].flavor = 'arc' $
                       else if imtype eq 'FLAT' then planstr[i].flavor = 'domeflat' $
                       else if imtype eq 'SKY' then planstr[i].flavor = 'twiflat' $
                       else if imtype eq 'STD' then planstr[i].flavor = 'std' $
                       else if imtype eq 'OBJECT' then planstr[i].flavor = 'science'

                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                       planstr[i].maskname = strtrim(sxpar(hdr, 'MASKNAME'))
                  ; Note that GRISNAME will only be defined for blue CCDs,
                  ; and GRANAME will only be defined for red CCDs. 
                  ; MSWAVE is only defined for red CCDs
                       instrument = strtrim(sxpar(hdr, 'INSTRUME'))
                       planstr[i].instrument = instrument
                       planstr[i].grating = strtrim(sxpar(hdr, 'GRATNAME'))
                       
                    endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*bok*') $
                    THEN BEGIN
                  ;;---------------------------
                  ;; This is the B&C spectrograph at the Bok 2.3-meter
                  ;; telescope; added by Moustakas on 2011-Jun-08
                  ;;---------------------------
                  exptime = strcompress(sxpar(hdr, 'EXPTIME'), /rem)
                  obstype = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                  object =  strcompress(sxpar(hdr, 'OBJECT'), /rem)
                  ccd    =  strcompress(sxpar(hdr, 'DETECTOR'), /rem)
                  instrument = 'Bok B&C'
                  mjd[i] = double(sxpar(hdr, 'MJD-OBS'))
                  dispers = sxpar(hdr, 'DISPERSE')
                  planstr[i].grating = strcompress(dispers, /rem) 
                  planstr[i].wave = strcompress(sxpar(hdr, 'TILTPOS'), /rem)
                  if obstype EQ 'flat' then planstr[i].flavor = 'domeflat' $
                  else if obstype EQ 'comp' then planstr[i].flavor = 'arc' $
                  else if obstype EQ 'zero'then planstr[i].flavor = 'bias' $
                  else if obstype EQ 'object' then begin
                     planstr[i].flavor = 'science' 
;                    if exptime LT 200.0 THEN planstr[i].flavor = 'std' $
;                    else planstr[i].flavor = 'science' 
                  endif else if obstype EQ 'focus' THEN $
                     planstr[i].flavor = 'focus' $
                  else message,  'Unrecognized obstype'
                  planstr[i].exptime = exptime
                  planstr[i].target = object
                  planstr[i].maskname =  $
                     strcompress(sxpar(hdr, 'APERTURE'), /rem)
                  planstr[i].instrument = instrument
              endif else if strmatch(sxpar(hdr, 'TELESCOP'), '*SOAR*') $
                then begin
                       ;;---------------------------
                       ;; This is the SOAR 4.1m
                       ;;---------------------------
                       exptime = sxpar(hdr, 'EXPTIME')
                       ;; Decide if this is twilight: sun_angle is the
                       ;; angle of the sun *above* the horizon, so must
                       ;; be less than -12 if darker than 12-degree twi.
                       ;observatory, 'ca', obs_struct
                       ;tzone = obs_struct.tz
                       ;lng = 360.d0 - obs_struct.longitude
                       ;lat = obs_struct.latitude
                       ;jd = sxpar(hdr, 'JUL-DATE')
                       ;sunpos, jd, ra1, dec1 ; returns degrees
                       ;zenpos, jd, ra2, dec2 ; returns radians
                       ;ra2 = ra2 * DRADEG
                       ;dec2 = dec2 * DRADEG
                       ;sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
                       ;caha_hh = caha_hdr(hdr)
                       ;imagetyp = strcompress(sxpar(hdr, 'IMAGETYP'), /rem)
                       
                       ; Use filename to find out what type of obs
                       
                       tmp_filetyp = STRSPLIT(filenames[i],'.',/EXTRACT)
                       filetyp = tmp_filetyp[1]

                       IF filetyp EQ 'Bias' OR filetyp EQ 'bias' $
                       THEN planstr[i].flavor = 'bias' $
                       ELSE IF filetyp EQ 'arc' OR filetyp EQ 'Arc' THEN $
                          planstr[i].flavor = 'arc' $
                       ELSE IF filetyp EQ 'flat' OR filetyp EQ 'InternalFlat' THEN $
                          planstr[i].flavor = 'iflat' $
                       ELSE IF filetyp EQ 'domeflat' OR filetyp EQ 'DomeFlat' THEN $
                         planstr[i].flavor = 'domeflat' $
                       ELSE IF filetyp EQ 'skyflat' OR filetyp EQ 'Skyflat' THEN $
                         planstr[i].flavor = 'twiflat' $
                       ELSE BEGIN
                           IF exptime LT 350.0 THEN $
                             planstr[i].flavor = 'std' $
                           ELSE planstr[i].flavor = 'science'
                       ENDELSE 
                       ;ENDIF ELSE message,  'Unrecognized obstype'
                       planstr[i].exptime = exptime
                       planstr[i].target = strtrim(sxpar(hdr, 'OBJECT'))
                       planstr[i].maskname = strtrim(sxpar(hdr, 'SLIT'))

                       instr1 = strcompress(sxpar(hdr, 'INSTRUME'), /rem)
                       ;path = strcompress(sxpar(caha_hh, 'PATH'), /rem)
                       instrument = instr1; + '-' + path
                       planstr[i].instrument = instrument
                       planstr[i].grating = strtrim(sxpar(hdr, 'GRATING'))
                       

                       
                                          
                   ENDIF ELSE $
                     splog, 'WARNING: Unknown instrument for ', filenames[i]
            ENDIF
            ENDFOR              ; End loop over files
           isort = sort(mjd)
           planstr = planstr[isort]
           logfile = repstr(planfile, '.par', '') + '.log'
           plotfile = repstr(planfile, '.par', '') + '.ps'
           
         hdr = ''
;         hdr = [hdr, '# Biases grouped by INSTRUMENT']
;         hdr = [hdr, '# Pixel flats grouped by INSTRUMENT+GRATING+MASKNAME+WAVE']
;         hdr = [hdr, '# Wavelength solutions from arcs grouped by INSTRUMENT+GRATING+MASKNAME+WAVE']
;         hdr = [hdr, '# Illumination flats grouped by INSTRUMENT+GRATING+MASKNAME+WAVE']
         hdr = [hdr, ' ']
         hdr = [hdr, "logfile '" + logfile + "'   # Log file"]
         hdr = [hdr, "plotfile '" + plotfile + "'   # Plot file"]
         hdr = [hdr, "indir '" + dirlist[idir] + "'   # Raw data directory"]
         hdr = [hdr, "tempdir Temp     # Temporary (working) directory"]
         hdr = [hdr, "scidir  Science  # Science output directory"]
;         if n_elements(maxobj) GT 0 then $
;            hdr = [hdr, "maxobj  "+string(maxobj,format='(i)') + $
;                        " # Maximum number of Objects per slit"] $
;         else 
         hdr = [hdr, "maxobj     5  # Maximum number of Objects per slit"] 
         hdr = [hdr, "minslit 20  # Minimum slit width"]
;         hdr = [hdr, "slitthresh 0.02 # Sets threshold for slit identification. 0.1 works best for Longslit spectra"]
         hdr = [hdr, "reduxthresh 0.01 # Sets the fraction of the brightest objects on each slit that is reduced"]
         hdr = [hdr, "idlutilsVersion '" + idlutils_version() $
                + "'  # Version of idlutils when building plan file"]
         hdr = [hdr, "LongslitVersion '" + longslit_version() $
                + "'  # Version of Longslit when building plan file"]
         
         ;; Write this plan file
         cd, olddir
         yanny_write, planfile, ptr_new(planstr), hdr = hdr, /align
         
      endif
   endfor ; End loop over directories

   cd, origdir

   return
end
;------------------------------------------------------------------------------
