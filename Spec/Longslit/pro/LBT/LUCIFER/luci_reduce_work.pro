PRO luci_reduce_work, afiles, bfiles, slitfile, scifile, waveqafile $
                      , FILESTD = FILESTD1 $
                      , TELLURIC = TELLURIC, WAVEFILE_TELL = WAVEFILE_TELL $
                      , CHK = CHK, WVCHK = WVCHK $
                      , pixflatfile = pixflatfile $
                      , illumflatfile = illumflatfile $
                      , darkfile = darkfile, SOFI = SOFI, CALIB = CALIB $
                      , NOHELIO = NOHELIO, PROF_NSIGMA = PROF_NSIGMA, SIZE_OBJMASK = SIZE_OBJMASK $
                      , reduxthresh = reduxthresh 


t0 = systime(1)
  

IF KEYWORD_SET(CHK) THEN set_plot, 'X'
;----------
; Set defaults
if (NOT keyword_set(box_rad)) then box_rad = 8L

slitmask = mrdfits(slitfile, 0)
tset_slits = mrdfits(slitfile,1)
IF KEYWORD_SET(FILESTD1) THEN FILESTD = FILESTD1
IF KEYWORD_SET(FILESTD) THEN BEGIN
    splog, 'Using telluric as crutch from ' + filestd
    stdstruct = xmrdfits(filestd, 3, /silent)
    stdmax = max(stdstruct.PEAKFLUX, stdind)
    stdtrace = stdstruct[stdind].XPOS
 ENDIF ELSE stdtrace = 0

diff = luci_skysub(afiles, bfiles, tset_slits, ivar = ivar $
                   , waveimg = waveimg, sky = sky $
                   , STDTRACE = STDTRACE $
                   , obj_pos = obj_pos, obj_neg = obj_neg $
                   , telluric = telluric, WAVEFILE_TELL = WAVEFILE_TELL $
                   , chk = chk, WVCHK = WVCHK $
                   , qafile = waveqafile $
                   , pixflatfile = pixflatfile $
                   , illumflatfile = illumflatfile, darkfile = darkfile $
                   , SIZE_OBJMASK = SIZE_OBJMASK $
                   , peakthresh = reduxthresh $
                   , SOFI = SOFI, CALIB = CALIB)
; Read in order set structure and create ordermask
IF KEYWORD_SET(SOFI) THEN plate_scale = 0.288D ELSE plate_scale = 0.250D
;----------
; Loop over objects and extract
IF KEYWORD_SET(OBJ_POS) THEN npos = n_elements(obj_pos) ELSE npos = 0
IF KEYWORD_SET(OBJ_NEG) THEN nneg = n_elements(obj_neg) ELSE nneg = 0
nobj = npos + nneg

final_struct = 0
;; Loop over positive objects and extract
IF npos GT 0 THEN BEGIN
   FOR iobj = 0L, npos-1L DO BEGIN
      extract_pos = gnirs_extract(diff, ivar, waveimg, (slitmask GT 0) $
                                  , sky, obj_pos[iobj], plate_scale $
                                  , SN_GAUSS = SN_GAUSS, TELLURIC = TELLURIC)
      final_struct = struct_append(final_struct, extract_pos)
   ENDFOR
ENDIF
IF nneg GT 0 THEN BEGIN
   FOR iobj = 0L, nneg-1L DO BEGIN
      extract_neg = gnirs_extract(-diff, ivar, waveimg, (slitmask GT 0) $
                                  , sky, obj_neg[iobj], plate_scale $
                                  , SN_GAUSS = SN_GAUSS, TELLURIC = TELLURIC)
      final_struct = struct_append(final_struct, extract_neg)
   ENDFOR
ENDIF

scihdr = xheadfits(afiles[0])
IF nobj GT 0 THEN BEGIN 
   isort = sort(final_struct.xfracpos)
   final_struct = final_struct[isort]
   final_struct.OBJID = lindgen(nobj) + 1L
   ;; compute heliocentric correction
   IF NOT KEYWORD_SET(NOHELIO) THEN long_helio, scihdr, final_struct
ENDIF
   

;----------
; Write output file
splog, 'Writing FITS file ', scifile
mwrfits, float(diff), scifile, scihdr, /create
mwrfits, float(ivar), scifile
mwrfits, float(waveimg), scifile
mwrfits, final_struct, scifile

IF nobj NE 0 THEN niri_plotsci, scifile $
                                , hard_ps = repstr(scifile, '.fits', '.ps') $
                                , box = keyword_set(TELLURIC), /LUCI

splog, 'Elapsed time = ', systime(1)-t0, ' sec'

RETURN
END


