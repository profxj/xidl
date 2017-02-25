;+
; NAME:
;   long_reduce_params
;
; PURPOSE:
;
;   Return parameters for longslit reduction. 
; CALLING SEQUENCE:
;  long_reduce_params,scihdr
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;  flg_skyfile:  Flag indicating whether multiple skyfiles exist
;                 (0=No, 1=Yes)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   11-Mar-2005  Written by JH
;-  
;-----------------------------------------------------------------------------
PRO long_reduce_params, scihdr, bin, skyfile = skyfile, anamorph = anamorph $
                        , bsp = bsp, SN_GAUSS = SN_GAUSS $
                        , SKYTRACE = SKYTRACE, SKYSAMPLE = SKYSAMPLE $
                        , NOSHIFT = NOSHIFT, NCCD = NCCD, FWHM = FWHM $
                        , PEAK_SMTH = PEAK_SMTH, FILESTD = FILESTD $
                        , BOX_RAD = BOX_RAD, FLG_SKYFILE = flg_skyfile

;;  What telescope are we using?
telescope = strcompress(sxpar(scihdr[*, 0], 'TELESCOP'), /rem)
instrument = strcompress(sxpar(scihdr[*, 0], 'INSTRUME'), /rem)
detector   =  strcompress(sxpar(scihdr[*, 0], 'DETECTOR'), /rem)
mask    = strtrim(sxpar(scihdr[*,0], 'SLITNAME'))
fpa     = strcompress(sxpar(scihdr[*, 0], 'FPA'), /rem)
dichroic = strcompress(sxpar(scihdr[*, 0], 'DICHNAME'), /rem)
flg_skyfile = 0


IF strmatch(instrument, 'LRIS*') THEN BEGIN
   box_rad = 7
    IF instrument EQ 'LRISBLUE' THEN BEGIN 
        grism = strtrim(sxpar(scihdr, 'GRISNAME'))
        mask = strtrim(sxpar(scihdr, 'SLITNAME'))
        CASE grism OF
            ;; used qso on reddest slit #5 mask 0127
            '1200/3400': skyfile = GETENV('LONGSLIT_DIR') + $
                                   '/calib/sky/lris_sky_blue_1200.sav'
            '600/4000': BEGIN
               IF strmatch(dichroic, '*560*') THEN $
                  skyfile = GETENV('LONGSLIT_DIR') + $
                            '/calib/sky/lris_sky_blue_600.sav' $
               ELSE IF strmatch(dichroic, '*460*') THEN $
                  skyfile = GETENV('LONGSLIT_DIR') + $
                            '/calib/sky/lris_sky_blue_600_4000_d460.sav' $
               ELSE IF strmatch(dichroic, '*500*') THEN $
                  skyfile = GETENV('LONGSLIT_DIR') + $
                            '/calib/sky/lris_sky_blue_600.sav' $
               ELSE skyfile = 0
                                ;skyfile = GETENV('LONGSLIT_DIR') + $
                                ;'/calib/sky/lris_sky_blue_600.sav'
            END
            '400/3400':  skyfile = GETENV('LONGSLIT_DIR') + $
                                   '/calib/sky/lris_sky_blue_400.sav'
            '300/5000':  skyfile = GETENV('LONGSLIT_DIR') + $
                                   '/calib/sky/lris_sky_blue_300-5000.sav'
            ELSE: message, 'ERROR: Unknown grism'
         ENDCASE
        ;stop
        anamorph = 1.0D
        bsp = 1.2/double(bin[1])
        NOSHIFT = 0
        ;; archived standards for longslit
        IF strmatch(mask, '*long*') THEN BEGIN
            CASE grism OF
                '1200/3400': BEGIN
                    filestd = GETENV('LONGSLIT_DIR') + $
                      '/calib/standards/std_trace/LRIS/std-blue-1200.fits.gz'
                    NOSHIFT = 1 
                    ;; no shift for 1200 line longslit. not enough sky cnts.
                END
                '600/4000': filestd = 0 
                '300/5000': filestd = 0
;filestd =  GETENV('LONGSLIT_DIR') + $
;                  '/calib/standards/std_trace/LRIS/LKPstd-blue-long-300-5000.fits.gz'
                ELSE: BEGIN
                    splog, 'ERROR: no archived standard for this setup'
                    filestd = 0
                END
            ENDCASE
        ENDIF
    ENDIF ELSE BEGIN
        ;; Some defaults for the red side
        SN_GAUSS = 4.0D ;; Minimum SNR to do fancy optimal extraction 
        grating_tilt = sxpar(scihdr, 'GRANGLE')
        grating = strtrim(sxpar(scihdr, 'GRANAME'))
        grangle = double(strtrim(sxpar(scihdr, 'GRANGLE')))
        mask    = strtrim(sxpar(scihdr, 'SLITNAME'))
        mswave = double(strtrim(sxpar(scihdr, 'MSWAVE')))
        IF NOT keyword_set(grating_tilt) then anamorph = 0.95 $
        ELSE anamorph = cos(!Pi/180 * (grating_tilt-59.)) / $
          cos(grating_tilt*!Pi/180.)
        bsp = 0.75D/double(bin[1]) 
        ;; changed to 0.75 (from 0.7) to try to fix big glitches in sky 
        ;; subtraction
        NOSHIFT = 0 ;; apply shifts to take out spatial flexure 
        ;; archived sky spectra to take out wavelength flexure
        ;; skytrace=1 (trace night sky lines to do better sky subtraction)
        ;; skytrace=0 (use traces of arc spectrum (since there are not enough 
        ;;             sky lines)
        CASE grating OF
            ;; used qso on reddest slit #5 mask 0127
            '300/5000': BEGIN
                skyfile = GETENV('LONGSLIT_DIR') + $
                  '/calib/sky/lris_sky_red_300_grangle_18.86_multi.sav'
                if n_elements(SKYTRACE) EQ 0 then skytrace = 1
            END
            '400/8500': BEGIN 
                skyfile = GETENV('LONGSLIT_DIR') + $
                  '/calib/sky/lris_sky_red_400.sav'
                if n_elements(SKYTRACE) EQ 0 then skytrace = 1
            END
            '600/5000': BEGIN 
                skyfile = GETENV('LONGSLIT_DIR') + $
                  '/calib/sky/lris_sky_red_600_5000.sav'
                if n_elements(SKYTRACE) EQ 0 then skytrace = 1
            END
            '600/7500': BEGIN 
               flg_skyfile = 1
                skyfile = GETENV('LONGSLIT_DIR') + $  ;; This is a root not a file
                  '/calib/sky/lris_sky_red_600_7500_waverange'
;                skyfile = GETENV('LONGSLIT_DIR') + $
;                  '/calib/sky/lris_sky_red_600_7500_grangle_29.54_long.sav'
                if n_elements(SKYTRACE) EQ 0 then skytrace = 1
            END
            '831/8200': BEGIN 
                skyfile = GETENV('LONGSLIT_DIR') + $
                  '/calib/sky/lris_sky_red_600_7500_grangle_29.54_long.sav'
                if n_elements(SKYTRACE) EQ 0 then skytrace = 1
;                stop ;; Need to put in a real file here for real data
            END
            ;; used longslit data for 2146-0754
            '600/10000': BEGIN
                if strmid(sxpar(scihdr,'DATE'),10) LT '2009-07-01' then begin
                    skyfile = GETENV('LONGSLIT_DIR') + $
                      '/calib/sky/lris_sky_red_600_10000_grangle_30.34_long.sav'
                    if n_elements(SKYTRACE) EQ 0 then skytrace = 1
                endif else skytrace = 0
            END
            '900/5500': BEGIN
                skyfile = GETENV('LONGSLIT_DIR') + $
                  '/calib/sky/lris_sky_red_900.sav'
                if n_elements(SKYTRACE) EQ 0 then skytrace = 1
            END
            '1200/7500': BEGIN 
                ;; For Jan_2007 LRIS data longslit has no sky lines so 
                ;; uncomment this line to disable flexure correction 
                ;; for the longslit
                ;;IF strmatch(mask, '*long*') THEN skyfile = 0 $
                ;;ELSE 
                skyfile = GETENV('LONGSLIT_DIR') + $ 
                  '/calib/sky/lris_sky_red_1200_grangle_34.82_multi.sav'
                if n_elements(SKYTRACE) EQ 0 then skytrace = 0 ;; not enough lines for this setup
            END
            '1200/9000': BEGIN 
                ;; For Jan_2007 LRIS data longslit has no sky lines so 
                ;; uncomment this line to disable flexure correction 
                ;; for the longslit
                ;;IF strmatch(mask, '*long*') THEN skyfile = 0 $
                ;;ELSE 
                skyfile = GETENV('LONGSLIT_DIR') + $ 
                  '/calib/sky/lris_sky_red_1200_grangle_34.82_multi.sav'
                if n_elements(SKYTRACE) EQ 0 then skytrace = 0 ;; not enough lines for this setup
            END
            ELSE: message, 'ERROR: Unknown grating'
        ENDCASE
        ;;use archived standards for longslit data
        IF strmatch(mask, '*long*') THEN BEGIN
            CASE grating OF
                ;; used qso on reddest slit #5 mask 0127
                '300/5000': filestd = GETENV('LONGSLIT_DIR') + $
       '/calib/standards/std_trace/LRIS/std-red-300_grangle_18.86.fits.gz'
                '600/10000': begin
                   if strmid(sxpar(scihdr, 'DATE'), 10) LT '2009-07-01' $
                   then begin
                      filestd = GETENV('LONGSLIT_DIR') + $
       '/calib/standards/std_trace/LRIS/std-red-600-10000_grangle_30.34.fits.gz'
                   endif else filestd = 0
                end
                ;; for now use 10000 trace for 7500 also. 
                '600/7500': begin
                   IF strmid(sxpar(scihdr, 'DATE'), 10) LT '2009-07-01' $
                   then begin
                      filestd = GETENV('LONGSLIT_DIR') + $
     '/calib/standards/std_trace/LRIS/std-red-600-10000_grangle_30.34.fits.gz' 
                   ENDIF ELSE filestd = 0
                end
;                '1200/7500': filestd = GETENV('LONGSLIT_DIR') + $
;                  '/calib/standards/std_trace/LRIS/std-red-1200_grangle_34.82.fits.gz'
                ELSE: BEGIN
                    splog, 'ERROR: no archived standard for this setup'
                    filestd = 0
                END
            ENDCASE
        ENDIF
    ENDELSE
ENDIF ELSE IF strcmp(telescope, 'Gemini-North') THEN BEGIN
    box_rad = 7
    nccd = 3L
    grating = sxpar(scihdr[*, 0], 'GRATING')
    grating_tilt = sxpar(scihdr[*, 0], 'GRTILT')
    if NOT keyword_set(grating_tilt) then anamorph = 1.45 $
    else anamorph = cos(!Pi/180 * (grating_tilt-59.))/cos(grating_tilt*!Pi/180.)
    NOSHIFT = 1
    wave = float(strcompress(sxpar(scihdr[*, 0], 'GRWLEN'), /rem))
    stddir = getenv('LONGSLIT_DIR') + '/calib/standards/std_trace/GMOS/'
    files = findfile(stddir + '*.fits*', count = ct)
    temp1 = stregex(files, 'std-...\.fits.gz', /extr)
    wave_std = double(strmid(temp1, 4, 3))
    diff = min(abs(wave-wave_std), jstd)
    filestd = files[jstd]
    ;; Skyfile
    case strmid(grating,0,4) of
        'R150': skyfile = GETENV('LONGSLIT_DIR') + $
                          '/calib/sky/gmosn_sky_R150.sav'
        'R400': skyfile = GETENV('LONGSLIT_DIR') + $
                          '/calib/sky/gmosn_sky_R400.sav'
        'R831': skyfile = GETENV('LONGSLIT_DIR') + $
                          '/calib/sky/gmosn_sky_R831.sav'
        'B600': skyfile = GETENV('LONGSLIT_DIR') + $
                          '/calib/sky/gmoss_sky_B600.sav'
        ELSE: skyfile=0
    endcase
ENDIF ELSE IF strcmp(telescope, 'Gemini-South') THEN BEGIN
    box_rad = 7
    nccd = 3L
    grating = sxpar(scihdr[*, 0], 'GRATING')
    grating_tilt = sxpar(scihdr[*, 0], 'GRTILT')
    if NOT keyword_set(grating_tilt) then anamorph = 1.45 $
    else anamorph = cos(!Pi/180 * (grating_tilt-59.))/cos(grating_tilt*!Pi/180.)
    NOSHIFT = 1
    wave = float(strcompress(sxpar(scihdr[*, 0], 'GRWLEN'), /rem))
    stddir = getenv('LONGSLIT_DIR') + '/calib/standards/std_trace/GMOS/'
    files = findfile(stddir + '*.fits*', count = ct)
    temp1 = stregex(files, 'std-...\.fits.gz', /extr)
    wave_std = double(strmid(temp1, 4, 3))
    diff = min(abs(wave-wave_std), jstd)
    filestd = files[jstd]
    ;; Skyfile
    case strmid(grating,0,4) of
        'B600': skyfile = GETENV('LONGSLIT_DIR') + $
                          '/calib/sky/gmoss_sky_B600.sav'
        'R400': skyfile = GETENV('LONGSLIT_DIR') + $
                          '/calib/sky/gmosn_sky_R400.sav' ;; This should be replaced [JXP]
        ELSE: skyfile=0
    endcase
ENDIF ELSE IF strmatch(telescope, '*DuPont*') THEN BEGIN
   box_rad = 7
   nccd = 1
   anamorph = 1
   bsp = 0.72
   skyfile = GETENV('LONGSLIT_DIR') + $
             '/calib/sky/dupont_sky_bcs_600.sav'
ENDIF ELSE IF strmatch(telescope, 'mmt') THEN BEGIN
    bsp = 0.72
    box_rad = 7
    PEAK_SMTH = 10.0D
    IF strmatch(detector, 'mmtblue*') then BEGIN
       grating = strcompress(sxpar(scihdr, 'DISPERSE'), /rem)
       IF strmid(sxpar(scihdr,'DATE'),10) LT '2014-10-01' THEN BEGIN
          CASE grating OF
             '300GPM':  skyfile = GETENV('LONGSLIT_DIR') + $
                                  '/calib/sky/mmt_sky_bcs_300.sav'
             '800GPM':  skyfile = GETENV('LONGSLIT_DIR') + $
                                  '/calib/sky/mmt_sky_bcs_800.sav'
             '1200GPM': skyfile = 0
             '500GPM': skyfile = 0
             ELSE: message, 'ERROR: Unknown grating'
          ENDCASE
       ENDIF ELSE BEGIN
          ;; KHRR -- for new detector, 23 Jan 2015
          CASE grating OF
             '300GPM':  skyfile = GETENV('LONGSLIT_DIR') + $
                                  '/calib/sky/mmt_sky_bcs_300_2015jan23.sav'
             '800GPM':  skyfile = GETENV('LONGSLIT_DIR') + $
                                  '/calib/sky/mmt_sky_bcs_800_2015jan23.sav'
             '1200GPM': skyfile = 0
             '500GPM': skyfile = GETENV('LONGSLIT_DIR') + $
                                  '/calib/sky/mmt_sky_bcs_500_2015jan23.sav'
             ELSE: message, 'ERROR: Unknown grating'
          ENDCASE
       ENDELSE 
    ENDIF ELSE BEGIN 
       grating = strcompress(sxpar(scihdr, 'DISPERSE'), /rem)
       CASE grating OF
          '600-6310':  skyfile = GETENV('LONGSLIT_DIR') + $
                                 '/calib/sky/mmt_sky_rcs_600.sav'
          ELSE: message, 'ERROR: Unknown grating'
       endcase
;       skyfile = 0
;       bsp = 0.72
    ENDELSE
    nccd = 1
    anamorph = 1
    NOSHIFT = 1
ENDIF ELSE IF strmatch(fpa, 'DBSP_BLUE') THEN BEGIN
   box_rad = 7
   nccd = 1
   anamorph = 1
   NOSHIFT = 1
   bsp = 0.72
ENDIF ELSE IF strmatch(instrument, 'DIS') THEN BEGIN
   box_rad = 7
   nccd = 1
   anamorph = 1
   NOSHIFT = 1
   IF strmatch(detector, '*red*') THEN bsp = 1.0d ELSE bsp = 0.7d
   IF strmatch(detector, '*red*') and n_elements(skytrace) EQ 0 $
   THEN skytrace = 1 ELSE skytrace = 0
ENDIF ELSE IF strmatch(instrument, '*IMACS*') THEN BEGIN
   box_rad = 7
   nccd = 1
   anamorph = 1
   NOSHIFT = 1
   bsp = 0.7D
   if n_elements(skytrace) EQ 0 then skytrace = 1
ENDIF ELSE IF (stregex(instrument,'.*kast.*',/boolean,/fold_case) eq 1) OR $
  (stregex(sxpar(scihdr,'VERSION'),'kast*', /bool, /fold_case) EQ 1) then begin
   box_rad = 7
   nccd = 1
   anamorph = 1
   NOSHIFT = 1
   bsp = 0.7D
   skytrace = 0
   if strtrim(sxpar(scihdr,'SPSIDE'),2) EQ 'blue' OR $
     strmid(sxpar(scihdr,'VERSION'),0,5) EQ 'kastb' then begin
       grating = strtrim(sxpar(scihdr,'GRISM_N'),2)
       CASE grating OF
           '452/3306':  skyfile = GETENV('LONGSLIT_DIR') + $
             '/calib/sky/kast_sky_blue_452.sav'
           '600/4310':  skyfile = GETENV('LONGSLIT_DIR') + $
             '/calib/sky/kast_sky_blue_600.sav'
           '830/3460':  skyfile = GETENV('LONGSLIT_DIR') + $
             '/calib/sky/kast_sky_blue_830.sav'
           ELSE: skyfile = 0
        ENDCASE
    endif else begin
        grating = strtrim(sxpar(scihdr,'GRATNG_N'),2)
        CASE grating OF
            '600/7500':  skyfile = GETENV('LONGSLIT_DIR') + $
             '/calib/sky/kast_sky_red_600_7500.sav'
            '830/8460':  skyfile = GETENV('LONGSLIT_DIR') + $
             '/calib/sky/kast_sky_red_830_8460.sav'
            '1200/5000':  skyfile = GETENV('LONGSLIT_DIR') + $
             '/calib/sky/kast_sky_red_1200_5000.sav'
            ELSE: skyfile = 0
         ENDCASE
    endelse
ENDIF ELSE IF  strmatch(telescope, '*kp4m*') THEN BEGIN
   box_rad = 3.5D
   FWHM = 1.8D ;; large pixels = 0.68"
   nccd = 1
   anamorph = 1.0
   NOSHIFT = 1
   bsp = 0.7D
   grating = strcompress(sxpar(scihdr, 'DISPERSE'), /rem)
   if grating eq 0 then grating = 'KPC10A'
   CASE grating OF
      'KPC10A':  skyfile = GETENV('LONGSLIT_DIR') + $
                           '/calib/sky/kpno_sky_rcspec_kpc-10a.sav'
      ELSE: message, 'ERROR: Unknown grating'
   ENDCASE
   if n_elements(skytrace) EQ 0 then skytrace = 0
ENDIF ELSE IF strmatch(telescope, '*bok*') THEN BEGIN ; jm11jun08ucsd
   box_rad = 9D ; 15" for 1.66666"/pixel                                                                                         
   FWHM = 2.0D
   peak_smth = 3D
   nccd = 1
   anamorph = 1.0
   NOSHIFT = 1
   bsp = 0.7D
;  if n_elements(skytrace) EQ 0 then skytrace = 0
   skyfile = getenv('LONGSLIT_DIR')+'/calib/sky/bok_sky_bcs_400.sav'
ENDIF ELSE IF  strmatch(telescope, '*CA-3.5*') THEN BEGIN
   caha_hh = caha_hdr(scihdr)
   path = strcompress(sxpar(caha_hh, 'PATH'), /rem)
   grating =  (path EQ 'BLUE' ? sxpar(caha_hh, 'GRAT1_NA') : '') $
              + (path EQ 'RED'  ? sxpar(caha_hh, 'GRAT2_NA') : '')
   CASE path OF
      'BLUE':  skyfile = GETENV('LONGSLIT_DIR') + $
                           '/calib/sky/CAHA_T13_BS6800_sky.sav'
      ELSE: message, 'ERROR: Unknown path'
   END
   box_rad = 7
   nccd = 1
   anamorph = 1
   NOSHIFT = 1
   bsp = 0.7D
   if n_elements(skytrace) EQ 0 then skytrace = 0
ENDIF ELSE IF strmatch(instrument, 'FORS2*') THEN BEGIN
   box_rad = 7
   anamorph = 1.0D
   grism = esopar(scihdr[*, 0], 'HIERARCH ESO INS GRIS1 NAME')
   IF strmatch(grism, '*GRIS_600B*') THEN BEGIN
      bsp = 1.0
      skysample = 0
   ENDIF ELSE IF strmatch(grism, '*GRIS_1200B*') THEN BEGIN
      skyfile = GETENV('LONGSLIT_DIR') + $
                '/calib/sky/fors2_sky_1200B.sav'
      bsp = 0.7
      skysample = 1  
      ;stop
   ENDIF ELSE IF strmatch(grism, '*GRIS_600RI*') THEN BEGIN
      skyfile = GETENV('LONGSLIT_DIR') + $
                '/calib/sky/fors2_sky_600RI.sav'
      bsp = 0.7
      skysample = 1  
   ENDIF ELSE IF strmatch(grism, '*GRIS_600V*') THEN BEGIN
      skyfile = GETENV('LONGSLIT_DIR') + $
                '/calib/sky/fors2_sky_600V.sav'
      bsp = 0.7
      skysample = 1 
   ENDIF ELSE BEGIN
      bsp = 0.7
      skysample = 1  
      ;; if skysample is set, determine optimal break-point spacing 
      ;; directly measuring how well we are sampling of the sky. The
      ;; bsp in this case correspons to the minimum distance between 
      ;; breakpoints which we allow. 
   ENDELSE 
   ;; FORS has not big flexure issue (apparently). Efforts to measure
   ;; it have failed. 
   NOSHIFT = 1 
   SN_GAUSS = 3.0D ;; Minimum SNR to do fancy optimal extraction 
   SKYTRACE = 0    ;; 
   filestd = 0
   nccd = 1
ENDIF ELSE IF strmatch(telescope, '*LBT-SX*') THEN BEGIN
   grating = strtrim(sxpar(scihdr, 'GRATINFO'),2)
   grating = strmid(grating,0,7)
   CASE grating OF 
       '250l/mm': skyfile = GETENV('LONGSLIT_DIR') + $
         '/calib/sky/mods_sky_red_670_long.sav' 
       '400l/mm': skyfile = GETENV('LONGSLIT_DIR') + $
         '/calib/sky/mods_sky_blue_450_long.sav'
       ELSE: skyfile = 0 ;message, 'ERROR: Unknown path'
   END 
   box_rad = 20
   anamorph = 1.0D
   bsp = 0.7
   skysample = 1
   fwhm = 12.0  ;; MODS has small pixels (0.12"/pix)
   ;; if skysample is set, determine optimal break-point spacing 
   ;; directly measuring how well we are sampling of the sky. The
   ;; bsp in this case correspons to the minimum distance between 
   ;; breakpoints which we allow. 
   NOSHIFT = 0
   ;;SN_GAUSS = 10.0D ;; Minimum SNR to do fancy optimal extraction 
   SKYTRACE = 0    ;; 
   filestd = 0
   nccd = 1
ENDIF ELSE BEGIN  ;; Unknown
   box_rad = 7
   nccd = 1
   anamorph = 1
   NOSHIFT = 1
   bsp = 0.7D
   if n_elements(skytrace) EQ 0 then skytrace = 0
ENDELSE
IF KEYWORD_SET(VERBOSE) THEN splog, 'Anamorphic factor is ', anamorph


RETURN
END 
