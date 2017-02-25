;+
; NAME:
;   tspec_waveimg
;
; PURPOSE:
;   Create plan file(s) for running the GNIRS pipeline
;
; CALLING SEQUENCE:
;   tspec_waveimg, [ indir, fileexpr, planfile= ]
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
;   tspecs_plan_struct()
;
; REVISION HISTORY:
;   Jan-2012  Written by JXP; Modified gnris_waveimg by JFH
;-
FUNCTION TSPEC_WAVEIMG, image, hdr, tset_slits $
                        , piximg = piximg, MXSHIFT = MXSHIFT $
                        , LINLIST = LINLIST, SIGREJ = SIGREJ $
                        , BOX_RAD = BOX_RAD, QAFILE = QAFILE, CHK = CHK
; Gain of GNIRS Aladdin II detector
;GAIN = 13.0  ; e-/ADU
if not keyword_set(MXSHFT) then MXSHFT = 120L
if not keyword_set(SIGREJ) then sigrej = 2.0D
if (NOT keyword_set(box_rad)) then box_rad = 8L
IF NOT KEYWORD_SET(SIG2D) THEN SIG2D = 3.5D

if not keyword_set(QAFILE) then qafile = 'test_wave.ps'

if not keyword_set(LINLIST) then $
  linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/TRIPLESPEC_OH.lst'

dims = size(image, /dimens)
nx = dims[0]
ny = dims[1]
dimt = size(tset_slits.coeff, /dimen)
norders = dimt[1]
slit_hdr = strtrim(sxpar(hdr, 'SLIT'))
slit_str = strsplit(slit_hdr, 'arcsec', /extr)
slit = double(slit_str)
slit = 1.
;stop
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
        left  = floor(trace[j, iorder] - BOX_RAD) > 0
        right = ceil(trace[j, iorder] + BOX_RAD) > BOX_RAD
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
order_vec = [7, 6, 5, 4, 3L] 


restore, file = getenv('XIDL_DIR') + $
         '/Spec/Longslit/calib/linelists/tspec_wave_archive.sav'
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
    ;; Generate a better continuum
    bkspace = 20.0
    rmed = median(arc1d[*,iorder], round(bkspace))
    ;x_splot, arc1d[*,iorder], ytwo=rmed, /bloc
    gdr = where(rmed NE 0.)
    ;x_splot, arc1d[gdr,iorder]/rmed[gdr], /bloc
    djs_iterstat, arc1d[gdr,iorder]/rmed[gdr], median=dmed, sigma=dsig 
    gda = where(arc1d[*,iorder] LT rmed+3*dsig)
    ;x_splot, arc1d[*,iorder], xtwo=gda, ytwo=arc1d[gda,iorder], /bloc
    flg = lonarr(ny)
    flg[gda] = 1
    spec_set = bspline_iterfit(findgen(ny), arc1d[*,iorder] $
                               , invvar = (arc1d[*,iorder] LE 1d6 and flg) $
                               , bkspace = bkspace $
                               , yfit = autofit $
                               , upper = upper, lower = lower, nord = 3 $
                               , maxrej = 10, outmask = outmask $
                               , /silent, /sticky)
    ;x_splot, findgen(ny), arc1d[*,iorder], ytwo=autofit, /block
;                               , invvar = (arc1d[*,iorder] GE (0.0 < median(arc1d[*,iorder]) $ ;; JXP
;                                                               AND arc1d[*,iorder] LE 1d6) $
    ;; Call
    x_templarc, arc1d[*, iorder], lines, calib[iorder], MSK = msk $
                , SHFT = shft, ALL_PK = all_pk, PKWDTH = pkwdth $
                , TOLER = TOLER, FORDR = 9, PKSIG = fsig1, FLG = flg_templ $
                , BEST = best, /THIN, AUTOFIT=autofit
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
    fsig_order = [5.0, 5.0, 5.0, 5.0, 4.0, 4.0] ; order 6 is more liberal
    ;fsig_order = [10.0, 10.0, 10.0, 10.0, 8.0, 4.0] 
;;  order 8 is more liberal
    fsig = fsig_order[iorder]   ; more refined line identification
    x_templarc, arc1d[*, iorder], lines, tmp_fit, MSK = msk $
                , PKWDTH = pkwdth, TOLER = TOLER, FORDR = 9, SHFT = shft $
                , PKSIG = fsig, /THIN, AUTOFIT=autofit ;, /FWEIGHT
;;  Check the number of good lines
    gdfit = where(lines.flg_plt EQ 1, ngd)
    fin_fit = calib[iorder]
    fit = x_fitrej(lines[gdfit].pix, lines[gdfit].wave, FITSTR = fin_fit $
                   , REJPT = rejpt, GDPT = gdpt, RMS=rms)
    print, order_vec[iorder], ': RMS (microns) = ', rms
    ngd =  n_elements(gdpt)
    sv_lines[iorder].nlin =  n_elements(gdpt)
    sv_lines[iorder].pix[0:ngd-1] =  lines[gdfit[gdpt]].pix
    sv_lines[iorder].wv[0:ngd-1] =  lines[gdfit[gdpt]].wave
  ;stop
    ;; Check
    ;CHK=1
    ;if keyword_set(CHK) and iorder EQ 4 then begin
    ;   npix = n_elements(arc1d[*,iorder])
    ;   restore, file = getenv('XIDL_DIR') + $
    ;            '/Spec/Longslit/calib/linelists/tspec_wave_archive.sav'
    ;   x_splot, x_calcfit(findgen(npix), fitstr=fin_fit), $
    ;            arc1d[*,iorder], $
    ;            xtwo=x_calcfit(findgen((size(archive_arc,/dim))[0]), fitstr=calib[iorder]), $
    ;            ytwo=archive_arc[*,iorder], /bloc, $
    ;            xthr=lines[gdfit].wave, ythr=fltarr(ngd)
    ;endif
ENDFOR

save_file = repstr(qafile, '.ps', '.sav')
save, sv_lines, arc1d, file = save_file

ordr_str = replicate(create_struct('FLG_ANLY', 0L, 'ORDER', 0L), norders)
ordr_str.flg_anly = 1L
ordr_str.ORDER = order_vec
;; nycoeff is along pixel direction for each order. This should 
;; be 3 since the maximum 1-d fits for any of the orders is 3. 
nycoeff = 4L                    ; 4L    
;; nocoeff is the order direction. 5 seems to give better rms.
nocoeff = 4L                    ;5L
dum = x_fit2darc('dum', ordr_str, save_file, nycoeff = nycoeff $
                 , nocoeff = nocoeff, /CLOBBER $
                 , DEBUG = debug, CHKRES = chkres, out_str = wave_struct $
                 , QAFIL = QAFILE, YVAL = yval, ORIG = orig, SIGR = SIG2D $
                 , sz = dims, pixrms = pixrms)


;pixset = long_wavepix(image, tset_slits, FWHM = FWHM, pkwdth = pkwdth $
;                      , toler = toler, CHK = CHK, ONLY_SLITS = [1, 2, 3, 4, 5])
  omask = tspec_ordermask(tset_slits)
  omask_trim = omask
  imask = (image GT -20.0 AND image LT 1d5 AND omask_trim GT 0)
  wset = tspec_wavepix(imask*image, tset_slits, fwhm = fwhm $    
                       , sig_thresh = sig_thresh $
                       , nsig = NSIG $
                       , pkwdth = PKWDTH $
                       , TOLER = TOLER $
                       , med_err = med_err $
                       , BOX_RADIUS = box_radius $
                       , CHK = chk $
                       , verbose = verbose $
                       , _EXTRA = keys ) 

  ;; 
  piximg = tspec_wpix2image(wset, tset_slits, wave_struct=wave_struct, waveimg=waveimg)
  ;piximg = gnirs_wpix2image(pixset, tset_slits, wave_struct = wave_struct $
  ;                        , waveimg = waveimg)

  ;; Apply heliocentric correction
  observ = strtrim(sxpar(hdr, 'OBSERVAT'),2) ;; APO vs Palomar
  case observ of
     'APO': begin
        mjd = strmid(sxpar(hdr,'DATE-OBS'), 0, 10)
        time = strmid(sxpar(hdr,'DATE-OBS'), 12)
     end
     else: begin
        mjd = strmid(sxpar(hdr,'UTSHUT'), 0, 9)
        time = sxpar(hdr, 'TIME')
        observ='Palomar'
     end
  endcase

  ;; 
  ra = sxpar(hdr, 'RA')
  dec = sxpar(hdr, 'DEC')
  equinox = sxpar(hdr, 'EQUINOX')
  mjd = x_setjdate(mjd, time)

  x_radec, ra, dec, rad, decd

  ;; Call helio
  helio = (-1.)*x_keckhelio(rad, decd, equinox, jd = mjd, OBS = observ)
  print, 'long_reduce: heliocentric correction :', $
         helio, ' km/s', format = '(a,f8.3,a)'
  hel_corr = sqrt( (1.d + helio/299792.458d) / $
                   (1.d - helio/299792.458d) )
  waveimg = waveimg*hel_corr
  ;stop

  RETURN, waveimg
END
