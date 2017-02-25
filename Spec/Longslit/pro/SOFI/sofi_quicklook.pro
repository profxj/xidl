path = '/Users/joe/DATA/SOFI_DATA/2011-09-22/'
TELLURIC = 0                  ;path + 'SOFI_0091.fits'
;;TELLURIC = path + 'SOFI_0099.fits'
;init = 107
;nseq = 1
init = 95
nseq = 4
IF KEYWORD_SET(TELLURIC) AND nseq NE 1 THEN message, 'Do not stack Tellurics'
even = 2*lindgen(nseq)
odd = 2*lindgen(nseq) + 1
anum = init + even
bnum = init + odd

afiles = path + 'SOFI_' + string(anum, FORMAT = '(I4.4)') + '.fits'
bfiles = path + 'SOFI_' + string(bnum, FORMAT = '(I4.4)') + '.fits'

IF KEYWORD_SET(TELLURIC) THEN scifile = 'tel-' + fileandpath(afiles[0]) $
ELSE scifile = 'sci-' + fileandpath(afiles[0]) 
;;wavefile = 'wave-' + fileandpath(afiles[0]) 
reduxpath = '//Users/joe/SOFI_pipeline/'
hdr = headfits(afiles[0])
object = strcompress(sxpar(hdr, 'OBJECT'), /rem)
targdir = reduxpath + '/' + object
IF FILE_TEST(targdir, /DIR) EQ 0 THEN spawn, 'mkdir ' + targdir
plate_scale = 0.288D
a_min_b = sofi_skysub(afiles, bfiles,ivar = ivar, sky_resids = sky_resids $
                      , hdr = scihdr, obj_pos = obj_pos, obj_neg = obj_neg $
                      , tset_slits = tset_slits, slitmask = slitmask $
                      , targdir = targdir $
                      , waveimg = waveimg, SKY_MODEL = SKY_MODEL $
                      , TELLURIC = TELLURIC, wavefile = wavefile $
                      , peakthresh = peakthresh, /CHK)
;; Kludged wavelengths
;ind_mode = WHERE(stregex(scihdr, 'HIERARCH ESO INS MODE', /bool))
;mode = strcompress(repstr(strmid(scihdr[ind_mode], 30, 15), "'", ''), /rem)
;CASE mode OF 
;   'LONG_SLIT_RED': BEGIN
;      ;; H+K
;      dlam = (2.52d4-1.53d4)/1024.0d
;      waveimg = 1.53d4 + waveimg*dlam
;   END
;   'LONG_SLIT_K': BEGIN
;      ;; K-band
;      dlam = (2.52d4-2.0d4)/1024.0d
;      waveimg = 2.0d4 + waveimg*dlam
;   END
;   'LONG_SLIT_H': BEGIN
;      ;; H-band
;      dlam = (1.80d4-1.50d4)/1024.0d
;      waveimg = 1.50d4 + waveimg*dlam
;   END
;   ELSE: message, 'Unrecognized setup'
;ENDCASE

; Loop over objects and extract
IF KEYWORD_SET(obj_pos) THEN npos = n_elements(obj_pos) ELSE npos = 0L
IF KEYWORD_SET(obj_neg) THEN nneg = n_elements(obj_neg) ELSE nneg = 0L
final_struct = 0
FOR iobj = 0L, npos -1L DO BEGIN
   extract = gnirs_extract(a_min_b-sky_resids, ivar, waveimg, slitmask $
                           , sky_model, obj_pos[iobj], plate_scale)
   final_struct = struct_append(final_struct, extract)
ENDFOR
FOR iobj = 0L, nneg -1L DO BEGIN
   extract = gnirs_extract(sky_resids-a_min_b, ivar, waveimg, slitmask $
                           , sky_model, obj_neg[iobj], plate_scale)
   final_struct = struct_append(final_struct, extract)
ENDFOR
isort = sort(final_struct.xfracpos)
final_struct = final_struct[isort]
final_struct.OBJID = lindgen(n_elements(final_struct)) + 1L

;----------
; Write output file
splog, 'Writing FITS file ', scifile
mwrfits, float(a_min_b), scifile, scihdr, /create
mwrfits, float(sky_resids), scifile
mwrfits, float(ivar), scifile
mwrfits, float(waveimg), scifile
mwrfits, final_struct, scifile

;IF npos NE 0 OR nneg NE 0 THEN $
;   niri_plotsci, scifile, hard_ps = repstr(scifile, '.fits', '.ps') $
;  , box = keyword_set(TELLURIC)

;splog, 'Elapsed time = ', systime(1)-t0, ' sec'


END
