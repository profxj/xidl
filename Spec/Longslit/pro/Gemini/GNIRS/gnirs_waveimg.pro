FUNCTION GNIRS_WAVEIMG, image, hdr, tset_slits $
                        , piximg = piximg, MXSHIFT = MXSHIFT $
                        , LINLIST = LINLIST, SIGREJ = SIGREJ $
                        , BOX_RAD = BOX_RAD, QAFILE = QAFILE, CHK = CHK
; Gain of GNIRS Aladdin II detector
GAIN = 13.0  ; e-/ADU
if not keyword_set(MXSHFT) then MXSHFT = 120L
if not keyword_set(SIGREJ) then sigrej = 2.0D
if (NOT keyword_set(box_rad)) then box_rad = 8L
IF NOT KEYWORD_SET(SIG2D) THEN SIG2D = 2.5D

if not keyword_set(LINLIST) then $
  linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/oh_linelist.lst' 

dims = size(image, /dimens)
nx = dims[0]
ny = dims[1]
dimt = size(tset_slits.coeff, /dimen)
norders = dimt[1]
slit_hdr = strtrim(sxpar(hdr, 'SLIT'))
slit_str = strsplit(slit_hdr, 'arcsec', /extr)
slit = double(slit_str[0])
pkwdth = 1.3*slit/0.15D         ; GNIRS plate scale is 0.15"
TOLER = pkwdth/3.0D
FWHM = pkwdth

; generate left and right edge of slits
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge

trace = (left_edge + right_edge)/2.0D
; Median filtering is more robust against cosmics
arc1d = fltarr(ny, norders)
FOR iorder = 0L, norders-1L DO BEGIN
    FOR j = 0L, ny-1L DO BEGIN
        left  = floor(trace[j, iorder] - BOX_RAD)
        right = ceil(trace[j, iorder] + BOX_RAD)
        sub_arc = image[left:right, j]
        djs_iterstat, sub_arc, median = median, sigrej = 2.0
        arc1d[j, iorder] = median
    ENDFOR
ENDFOR

;; Open line list
x_arclist, linlist, lines

fittmp = { fitstrct }
all_arcfit = replicate(fittmp, norders)
lintmp = { $
           pix: dblarr(90), $
           wv: dblarr(90), $
           nlin: 0 $
         }
sv_lines = replicate(lintmp, norders)
order_vec = [3, 4, 5, 6, 7, 8]


restore, file = getenv('XIDL_DIR') + $
         '/Spec/Longslit/calib/linelists/gnirs_wave_archive.sav'
FOR iorder = 0L, norders-1L DO BEGIN

;;  Zero out
    lines.flg_plt = 0
;;  Cross-correlate to determine shift. Take out average trend in arc
    step = lindgen(2*MXSHFT) - MXSHFT 
    arc_obj1 = djs_median(arc1d[*, iorder], width = 30, boundary = 'reflect')
    arc_ref1 = djs_median(archive_arc[*, iorder], width = 30 $
                          , boundary = 'reflect')
    sig_res = 4.0D
    nhalf =  long(sig_res)*4L
    xkern = dindgen(2*nhalf+1)-nhalf
    kernel = gauss1(xkern, [0.0, sig_res, 1.0])
    arc_obj = convol(arc_obj1, kernel, /edge_truncate)
    arc_ref = convol(arc_ref1, kernel, /edge_truncate)
    arc_obj_corr = (arc1d[*, iorder] - arc_obj) <  1.0d6
    arc_ref_corr = (archive_arc[*, iorder] - arc_ref) <  1.0d6
    corr = c2_correlate(arc_obj_corr, arc_ref_corr, step, /double)
    xcen = long_find_nminima(-corr, nfind = 1, width = 5 $
                        , minsep = 1, ypeak = ypeak1, npeak = npeak)
    xcen = double(xcen)      
    shft = -step[long(round(xcen))]
    print, iorder+3, shft, FORMAT = '(%"Order: %d Wave Shift: %d pixels")' 
;   Re-identify lines
    fsig_order1 = [7.0, 7.0, 7.0, 7.0, 7.0, 4.0]
;    fsig_order1 = [5.0, 5.0, 5.0, 5.0, 4.0, 4.0] ; order 6 is more liberal
    fsig1 = fsig_order1[iorder] ; more refined line identification
    x_templarc, arc1d[*, iorder], lines, calib[iorder], MSK = msk $
                , SHFT = shft, ALL_PK = all_pk, PKWDTH = pkwdth $
                , TOLER = TOLER, FORDR = 9, PKSIG = fsig1, FLG = flg_templ $
                , BEST = best, /THIN, /FWEIGHT
    if flg_templ EQ -1 then begin
        print, 'x_fitarc: Insufficient lines for AUTO!!'
        stop
    endif
    
;; Check the number of good lines
    gdfit = where(lines.flg_plt EQ 1, ngd)
    
;; FIRST (crude) FIT
    tmp_fit = {fitstrct}
    copy_struct, calib[iorder], tmp_fit, $
                 EXCEPT_TAGS = ['FFIT', 'NRM', 'RMS']
    tmp_fit.flg_rej = 1 
    tmp_fit.niter = 3 
    tmp_fit.maxrej = (10 >  ngd/2) 
    tmp_fit.minpt = 5
    IF iorder EQ 5 THEN tmp_fit.hsig = 3.0D $
    ELSE tmp_fit.hsig = sigrej
    IF iorder EQ 5 THEN tmp_fit.lsig = 3.0D $
    ELSE tmp_fit.lsig = sigrej
    tmp_fit.nord = calib[iorder].nord 
; (3 < calib[iorder].nord) ; Set first fit to 3
    fin_fit = tmp_fit
    fin_fit.nord = calib[iorder].nord
    fit = x_fitrej(lines[gdfit].pix, lines[gdfit].wave, FITSTR = tmp_fit $
                   , REJPT = rejpt, gdpt = gdpt)
    if fit[0] EQ -1 then begin
        print, 'x_fitarc: AUTO Failed!!'
        print, 'Retrying with 2nd order'
        tmp_fit.nord = 2
        fit = x_fitrej(lines[gdfit].pix, lines[gdfit].wave, FITSTR = tmp_fit $
                       , REJPT = rejpt, gdpt = gdpt)
    endif
    
;;      Grab new lines
    lines.flg_plt = 0
    shft = 0L                   ; No shift expected now!
;;      AutoID lines
;;    fsig_order = [5.0, 5.0, 5.0, 5.0, 4.0, 4.0] ; order 6 is more liberal
    fsig_order = [10.0, 10.0, 10.0, 10.0, 8.0, 4.0] 
;;  order 8 is more liberal
    fsig = fsig_order[iorder]   ; more refined line identification
    x_templarc, arc1d[*, iorder], lines, tmp_fit, MSK = msk $
                , PKWDTH = pkwdth, TOLER = TOLER, FORDR = 9, SHFT = shft $
                , PKSIG = fsig, /THIN, /FWEIGHT
;;  Check the number of good lines
    gdfit = where(lines.flg_plt EQ 1, ngd)
    fin_fit = calib[iorder]
    fit = x_fitrej(lines[gdfit].pix, lines[gdfit].wave, FITSTR = fin_fit $
                   , REJPT = rejpt, GDPT = gdpt)
    ngd =  n_elements(gdpt)
    sv_lines[iorder].nlin =  n_elements(gdpt)
    sv_lines[iorder].pix[0:ngd-1] =  lines[gdfit[gdpt]].pix
    sv_lines[iorder].wv[0:ngd-1] =  lines[gdfit[gdpt]].wave
ENDFOR

save_file = repstr(qafile, '.ps', '.sav')
save, sv_lines, file = save_file


ordr_str = replicate(create_struct('FLG_ANLY', 0L, 'ORDER', 0L), norders)
ordr_str.flg_anly = 1L
ordr_str.ORDER = order_vec
;; nycoeff is along pixel direction for each order. This should 
;; be 3 since the maximum 1-d fits for any of the orders is 3. 
nycoeff = 3L                    ; 4L    
;; nocoeff is the order direction. 5 seems to give better rms.
nocoeff = 5L                    ;5L
dum = x_fit2darc('dum', ordr_str, save_file, nycoeff = nycoeff $
                 , nocoeff = nocoeff, /CLOBBER $
                 , DEBUG = debug, CHKRES = chkres, out_str = wave_struct $
                 , QAFIL = QAFILE, YVAL = yval, ORIG = orig, SIGR = SIG2D $
                 , sz = dims, pixrms = pixrms)
IF pixrms[5] GE 1.0D THEN BEGIN
; This is a hack until I get better solutions for order 8
; order 8 wavelengths will just be junk if the pixel rms is too large. 
    sv_lines[5].nlin = 0
    save, sv_lines, file = save_file
    dum = x_fit2darc('dum', ordr_str, save_file, nycoeff = nycoeff $
                     , nocoeff = nocoeff, /CLOBBER $
                     , DEBUG = debug, CHKRES = chkres, out_str = wave_struct $
                     , QAFIL = QAFILE, YVAL = yval, ORIG = orig, SIGR = SIG2D $
                     , sz = dims, pixrms = pixrms)
ENDIF

;;chk=1
pixset = long_wavepix(image, tset_slits, FWHM = FWHM, pkwdth = pkwdth $
                      , toler = toler, CHK = chk $
                      , med_err = [0.1, 0.1, 0.1, 0.1, 0.1, 0.3])
;;, ISLIT = [1, 2, 3, 4, 5])
;; Trying to trace order 8 causes many problems so don't try for now. 
;; Could improve this by using the arcs. 
;pixset[5].COEFF2D = 0.0
;pixset[5].COEFF2D_INV = 0.0
piximg = gnirs_wpix2image(pixset, tset_slits, wave_struct = wave_struct $
                          , waveimg = waveimg)
; Apply heliocentric correction
ra = sxpar(hdr, 'RA')
dec = sxpar(hdr, 'DEC')
equinox = sxpar(hdr, 'EQUINOX')
jd = 2400000.5D + sxpar(hdr, 'MJD_OBS')
observatory = sxpar(hdr, 'OBSERVAT')
IF strmatch(observatory, '*South*') THEN BEGIN
   latitude =  -30.240750D
   longitude = 70.736693D
   altitude = 2722.0D
   helio = (-1.0d)*x_keckhelio(RA, DEC, EQUINOX, jd = JD $
                               , longitude = longitude, latitude = latitude $
                               , altitude = altitude)
ENDIF ELSE BEGIN 
   helio = (-1.0d)*x_keckhelio(RA, DEC, EQUINOX, jd = JD)
ENDELSE
hel_corr = sqrt( (1.d + helio/299792.458) / (1.d - helio/299792.458) )
print, 'gnirs_waveimg: Heliocentric correction applied -- ', helio, hel_corr
waveimg = waveimg*hel_corr

RETURN, waveimg
END
