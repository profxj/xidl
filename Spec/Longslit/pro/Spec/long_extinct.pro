;+
; NAME:
;  long_extinct
;  Version 1.1
;
; PURPOSE:
;  Return the extinction as a function of wavelength as a function of
;  observatory.  Also pass back the AIRMASS parsed from the header
;
; CALLING SEQUENCE:
;  LONG_EXTINCT, scihdr, WAVE_EXT = , MAG_EXT =
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;  WAVE_EXT=  -- Array of wavelengths where the extinction is
;                evaluated
;  MAG_EXT1=  -- Array of extinction values in magnitudes.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;  x_gsmooth
;
;
; REVISION HISTORY:
;-
FUNCTION LONG_EXTINCT, wave, scihdr, NOEXTINCT = NOEXTINCT $
                       , AIRMASS = AIRMASS, EXPTIME = EXPTIME1 $
                       , MAG_EXT = MAG_EXT1, NOTIME = notime

  compile_opt strictarr

  instrument = sxpar(scihdr, 'INSTRUME')
  longslit_dir = getenv('LONGSLIT_DIR')
; read in atmospheric extinction data
IF strmatch(sxpar(scihdr, 'TELESCOP'), 'Keck*') OR $
   strmatch(instrument, '*LRIS*') THEN BEGIN
    extinctfile = longslit_dir + '/calib/extinction/mkoextinct.dat'
    exptime = double(sxpar(scihdr, 'ELAPTIME'))
    airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE IF strmatch(sxpar(scihdr, 'OBSERVAT'), 'Gemini*') THEN BEGIN
   IF strmatch(sxpar(scihdr, 'TELESCOP'), 'Gemini-North') THEN $
      extinctfile = longslit_dir + '/calib/extinction/mkoextinct.dat' $
   ELSE IF strmatch(sxpar(scihdr, 'TELESCOP'), 'Gemini-South') THEN $
      extinctfile = longslit_dir + '/calib/extinction/ctioextinct.dat'
   exptime = double(sxpar(scihdr, 'ELAPSED')) > $
             double(sxpar(scihdr, 'EXPTIME'))
   airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE IF strmatch(sxpar(scihdr, 'DETECTOR'), 'mars*') OR $
   strmatch(sxpar(scihdr, 'DETECTOR'), '*ccd34*') OR $
   strmatch(sxpar(scihdr, 'TELESCOP'), 'mmt*') OR  $
   strmatch(sxpar(scihdr, 'TELESCOP'), '*kp4m*') OR $
   strmatch(sxpar(scihdr, 'TELESCOP'), '*3.5m*') OR $
   strmatch(sxpar(scihdr, 'TELESCOP'), '*bok*') OR $ ; jm11jun08ucsd
   strmatch(sxpar(scihdr, 'TELESCOP'), '*CA-3.5*') OR $
   strmatch(sxpar(scihdr, 'TELESCOP'), '*CA-2.2*') OR $
   strmatch(sxpar(scihdr, 'TELESCOP'), '*LBT-SX*') OR $
   strmatch(sxpar(scihdr, 'TELESCOP'), '*ESO-VLT-U1*') OR $
   strmatch(sxpar(scihdr, 'TELESCOP'), '*ESO-NTT*') $
THEN BEGIN
    extinctfile = longslit_dir + '/calib/extinction/kpnoextinct.dat'
    exptime = double(sxpar(scihdr, 'OBSTIME')) > double(sxpar(scihdr, 'EXPTIME'))
    airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE if (stregex(sxpar(scihdr,'INSTRUME'),'.*kast.*',$
                       /boolean,/fold_case) eq 1) or $
  (stregex(sxpar(scihdr,'VERSION'),'kast*', /bool, /fold_case) EQ 1) then begin
    extinctfile = longslit_dir + '/calib/extinction/mthamextinct.dat'
    exptime = double(sxpar(scihdr, 'EXPOSURE') > sxpar(scihdr,'EXPTIME'))
    ;; Calculate airmass
    ras = sxpar(scihdr, 'RA')
    decs = sxpar(scihdr, 'DEC')
    x_radec, ras, decs, rad, decd
    hangl = float(strsplit(sxpar(scihdr,'HA'),':',/extrac))
    if hangl[0] LT 0. then hangl = hangl[0]-hangl[1]/60. $
    else hangl = hangl[0]+hangl[1]/60.
    airmass = airmass(40., decd, hangl)
ENDIF ELSE IF  strmatch(strtrim(sxpar(scihdr, 'TELID'), 2), '200') THEN BEGIN
    print, 'Warning:  Using KPNO for Palomar!'
    extinctfile = longslit_dir + '/calib/extinction/kpnoextinct.dat'
    exptime = double(sxpar(scihdr, 'EXPTIME'))
    airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE IF strmatch(strtrim(sxpar(scihdr, 'INSTRUME'), 2), 'OSIRIS') THEN BEGIN
    extinctfile = longslit_dir + '/calib/extinction/lapalmaextinct.dat'
    exptime = double(sxpar(scihdr, 'EXPTIME'))
    airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE IF strmatch(sxpar(scihdr, 'INSTRUME'), '*ISIS*') THEN BEGIN
    extinctfile = longslit_dir + '/calib/extinction/lapalmaextinct.dat'
    exptime = double(sxpar(scihdr, 'EXPTIME'))
    airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE IF strmatch(sxpar(scihdr, 'INSTRUME'), '*EFOSC*') THEN BEGIN
    extinctfile = longslit_dir + '/calib/extinction/ctioextinct.dat'
    exptime = double(sxpar(scihdr, 'EXPTIME'))
    ind_am = WHERE(stregex(scihdr, 'HIERARCH ESO TEL AIRM END', /bool))
    airmass = double(strmid(scihdr[ind_am], 30, 14))
ENDIF ELSE IF strmatch(sxpar(scihdr, 'TELESCOP'), 'DuPont') THEN BEGIN
    extinctfile = longslit_dir + '/calib/extinction/ctioextinct.dat'
    exptime = double(sxpar(scihdr, 'EXPTIME'))
    airmass = double(sxpar(scihdr, 'AIRMASS'))
ENDIF ELSE message, 'Telescope not found'
IF NOT KEYWORD_SET(exptime) and not keyword_set(NOTIME) THEN $
  message, 'Could not find exposure time'

IF KEYWORD_SET(extinctfile) THEN begin
;   print, 'long_extinct: Using extinction file ', extinctfile
   readcol, extinctfile, wave_ext, mag_ext, format = 'F,F' , /silent
endif

EXPTIME1 = exptime

wave_min = min(wave_ext, jmin)
wave_max = max(wave_ext, jmax)
mag_ext1 = dblarr(n_elements(wave))
ext = dblarr(n_elements(wave)) + 1.0D
IF KEYWORD_SET(MAG_EXT) and not keyword_set(NOEXTINCT) THEN BEGIN
    inds = WHERE(wave GE wave_min AND wave LE wave_max)
    mag_ext1[inds] = interpol(mag_ext, wave_ext, wave[inds])
    linds = WHERE(wave LT wave_min, nl)
    IF nl GT 0 THEN mag_ext1[linds] = mag_ext[jmin]
    rinds = WHERE(wave GT wave_max, nr)
    IF nr GT 0 THEN mag_ext1[rinds] = mag_ext[jmax]
    if keyword_set(NOTIME) then exptime = 1.
    ext = 10.0D^(0.4D*mag_ext1*airmass)/exptime  ;; TAKE OUT EXPTIME TOO
 ENDIF
;IF NOT KEYWORD_SET(EXPTIME)  THEN exptime = 1.
;IF NOT KEYWORD_SET(AIRMASS)  THEN airmass = 1.

RETURN, EXT
END
