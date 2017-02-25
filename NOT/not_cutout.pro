PRO NOT_CUTOUT, datafile, image = image $
                , superbias = superbias, delta_bias = delta_bias $
                , superflat = superflat, superdark = superdark $
                , instbiasfile = instbiasfile $
                , psfile = psfile, bar = bar $
                , NAME = NAME, LOBAR = LOBAR $
                , LOWER = LOWER, UPPER = UPPER $
                , SCALE = SCALE, BRIGHTNESS = BRIGHTNESS, CONTRAST = CONTRAST $
                , hdr = hdr
   
IF n_params() LT 1 THEN BEGIN
  print, 'Usage: not_cutout,file'
ENDIF

common xatv_color, r_vector, g_vector, b_vector
common xatv_state, state, tot_state, nimg, flg_img

IF NOT KEYWORD_SET(SCALE) THEN SCALE = 1L
IF NOT KEYWORD_SET(EXTEN) THEN EXTEN = 0
ncut = SCALE*256L
plate_scale = 0.217D

IF NOT KEYWORD_SET(SUPERDARKFILE) THEN $
   superdarkfile = '/Users/jhennawi/NOT_run/redux/superdark.fits'
IF NOT KEYWORD_SET(SUPERBIASFILE) THEN $
   superbiasfile = '/Users/jhennawi/NOT_run/redux/superbias.fits' 
IF NOT KEYWORD_SET(SUPERFLATFILE) THEN $
   superflatfile = '/Users/jhennawi/NOT_run/redux/superflat_g.fits' 

;; Read in superbias image
IF NOT KEYWORD_SET(SUPERBIAS) THEN BEGIN
   IF NOT KEYWORD_SET(SUPERBIASFILE) THEN $
      message, 'Currently superbiasfile required'
   superbias = xmrdfits(superbiasfile, 0, superbiashdr $
                        , silent = (keyword_set(verbose) EQ 0))
   delta_bias = xmrdfits(superbiasfile, 1, silent = (keyword_set(verbose) EQ 0))
ENDIF

IF NOT KEYWORD_SET(SUPERDARK) AND KEYWORD_SET(superdarkfile) THEN $
   superdark = xmrdfits(superdarkfile, 0, superdarkhdr $
                        , silent = (keyword_set(verbose) EQ 0))
IF NOT KEYWORD_SET(SUPERFLAT) AND keyword_set(superflatfile) THEN $
   superflat = xmrdfits(superflatfile, 0, biashdr, $
                        silent = (keyword_set(verbose) EQ 0))

not_proc, datafile, image, ivar, superbias = superbias $
          , delta_bias = delta_bias, instbiasfile = instbiasfile $
          , superdark = superdark, mask = mask, superflat = superflat $
          , ihdr = hdr

dims = size(image, /dim)
nx = dims[0]
ny = dims[1]

IF NOT KEYWORD_SET(LOWER) THEN LOWER = 2.0
IF NOT KEYWORD_SET(UPPER) THEN UPPER = 8.0

skymode = fltarr(4)
skysig  = fltarr(4)
gaintweak = fltarr(4)
FOR chip = 1L, 4L DO BEGIN
   ipix = WHERE(mask EQ chip)
   sky, image[ipix], skymode1, skysig1, /SILENT
   skymode[chip-1] = skymode1
   skysig[chip-1]  = skysig1
ENDFOR
mean_sky = total(skymode)/4.0D
mean_sig = total(skysig)/4.0D
;; Force the images to have the same average sky value. 
;; This accounts for errors in the gain. 
FOR chip = 1L, 4L DO BEGIN
   ipix = WHERE(mask EQ chip)
   image[ipix] = image[ipix]*(mean_sky/skymode[chip-1])
ENDFOR

min = mean_sky - lower*mean_sig
max = mean_sky + upper*mean_sig

xatv, image, min = min, max = max ;, /LOG

x_cen = 1400
y_cen = 600
state.centerpix = [x_cen, y_cen]
IF state.invert_colormap EQ 0 THEN BEGIN
    state.invert_colormap = abs(state.invert_colormap - 1)
    r_vector = reverse(r_vector)
    g_vector = reverse(g_vector)
    b_vector = reverse(b_vector)
ENDIF

; Good for log
;state.brightness = 0.75
;state.contrast   = 0.20
; Good for linear
IF KEYWORD_SET(BRIGHTNESS) THEN state.brightness = brightness $
ELSE state.brightness = 0.30
IF KEYWORD_SET(CONTRAST) THEN state.contrast = contrast $
ELSE state.contrast = 0.15
xatv_stretchct, state.brightness, state.contrast
xatv_refresh


IF KEYWORD_SET(PSFILE) THEN xatv_joe_writeps, psfile, /CLOBBER

RETURN
END
