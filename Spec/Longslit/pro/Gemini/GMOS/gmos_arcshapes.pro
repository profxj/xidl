PRO GMOS_ARCSHAPES, gsfile, idfile, waveimg, slitfile $
                    , FWHMSET = FWHMSET $
                    , REGION = REGION, MAXSEP = MAXSEP, FWTOL = FWTOL $
                    , FWCOEFF = FWCOEFF, QAFILE = QAFILE

IF KEYWORD_SET(QAFILE) THEN dfpsplot, qafile, /color
IF NOT KEYWORD_SET(MAXSEP) THEN MAXSEP = 10L
IF NOT KEYWORD_SET(REGION) THEN REGION = 5L
IF NOT KEYWORD_SET(FWTOL)  THEN FWTOL = 0.5
IF NOT KEYWORD_SET(FWCOEFF) THEN FWCOEFF = 3L
if NOT keyword_set(box_rad) then box_rad = 5L


; parse the idfile
lines = djs_readlines(idfile)
str = strsplit(lines[6], 'features', /extract)
features = long(str[1])
; read in the identified lines
readcol, idfile, xpix_iraf, wave_img, wave, weight, ind1, mask, skipline = 7 $
         , numline = features, format = 'F,F,F,F,I,I'

;xpix = double(ny)-(xpix_iraf -1.0D)
good  = WHERE(mask NE 0, ngood)
wv_lines = wave[good]
fwhm = fltarr(ngood)
;xpix = xpix[good]

; read in the iraf mosaiced wavelength image
arcimg =  reverse(transpose(mrdfits(gsfile, 2)), 2)
;arcimg =  reverse(transpose(mrdfits(gsfile, 1)), 2)
dims = size(arcimg)
nx = dims[1]
ny = dims[2]

;----------
; Read in slit structure 
tset_slits = xmrdfits(slitfile, 1, silent = (keyword_set(verbose) EQ 0))

;; ------
;; Expand slit set to get left and right edge
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge

if (size(left_edge, /n_dimen) EQ 1) then nslit = 1 $
else nslit = (size(left_edge, /dimens))[1]
trace = (left_edge + right_edge)/2.0D

; compute 1-d arc spectrum
arc1d = fltarr(ny, nslit)
wav1d = fltarr(ny, nslit)
for islit = 0L, nslit-1L do begin
    splog, 'Working on slit #', islit+1, ' (', islit+1, ' of ', nslit, ')'
    ;; iteratively compute mean arc spectrum which is robust against cosmics
    FOR j = 0L, ny-1L DO BEGIN
        left  = floor(trace[j, islit] - BOX_RAD)
        right = ceil(trace[j, islit] + BOX_RAD)
        sub_arc  = arcimg[left:right, j]
        sub_wav = waveimg[left:right, j]
        djs_iterstat, sub_arc $
                      , mean = mean1, sigrej = 3.0, mask = mask1
        arc1d[j, islit] = mean1
        djs_iterstat, sub_wav $
                      , mean = mean2, sigrej = 3.0, mask = mask2
        wav1d[j, islit] = mean2
    ENDFOR
    xpix = interpol(lindgen(ny), wav1d[*, islit], wv_lines)
    ;; Compute and fit for the FWHM 
    FOR j = 0L, ngood-1L DO $
      fwhm[j] = long_arc_fwhm(arc1d[*, islit], xpix[j], maxsep)
    inmask = fwhm GT 0.0 
    ;; Now fit a polynomial to the fwhm and reject outliers
    fwhmset_temp1 = jfh_fit_reject(xpix, fwhm, fwcoeff, TOL = FWTOL $
                                   , outmask = fwhmmask, inmask = inmask $
                                   , FUNC = 'poly', YFIT = FWFIT $
                                   , XMIN = 0, XMAX = ny-1L)
    good_lines = WHERE(fwhmmask GT 0, ngline)
    IF ngline NE 0 THEN fwhm_med = djs_median(fwhm[good_lines]) $
    ELSE fwhm_med = 0.0
    fwhmset_temp = struct_addtags(fwhmset_temp1 $
                                  , create_struct('MEDIAN', fwhm_med))
    ;; --------------
    ;; Output some QA
    IF keyword_set(QAFILE) THEN BEGIN 
        ;; create some bogus inputs for the QA routine
        fit = {fitstrct}
        fit.flg_rej = 1 
        fit.niter = 5
        fit.maxrej = ny/2
        fit.minpt  = ny/2
        fit.hsig = 3.0
        fit.lsig = 3.0
        fit.nord = 10
        fit.FUNC = 'CHEBY'
        lines = replicate(create_struct('PIX', 0.0D, 'WAVE', 0.0D), ngood)
        lines.pix = xpix
        lines.wave = wv_lines
        wfit = x_fitrej(findgen(ny), wav1d[*, islit], FITSTR = fit $
                       , REJPT = rejpt)
        long_waveqa, lines, fit, arc1d[*, islit], -1, islit $
                     , fwhm, fwhmset_temp, fwhmmask, QAFILE
        !p.multi = 0
    ENDIF
    
    IF NOT KEYWORD_SET(fwhmset) THEN BEGIN
        fwhmset = replicate(fwhmset_temp, nslit)
        struct_assign, {junk:0}, fwhmset ; Zero-out all elements
    ENDIF
    copy_struct_inx, fwhmset_temp, fwhmset, index_to = islit
ENDFOR

IF KEYWORD_SET(QAFILE) THEN dfpsclose

RETURN
END
