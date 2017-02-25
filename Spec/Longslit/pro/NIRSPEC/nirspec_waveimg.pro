FUNCTION NIRSPEC_WAVEIMG, image, tset_slits,hdr $
                          , OUTFILE=OUTFILE $
                          , piximg = piximg, pixset=pixset $
                          , MXSHFT = MXSHFT, FILTER=filter $
                          , SIGREJ = SIGREJ $
                          , BOX_RAD = BOX_RAD, QAFILE = QAFILE, CHK = CHK

IF NOT KEYWORD_SET(PKSIG) THEN pksig = 4.0D
if not keyword_set(MXSHFT) then MXSHFT = 15L
if not keyword_set(FILTER) then filter = ''
if not keyword_set(SIGREJ) then sigrej = 2.
if (NOT keyword_set(box_rad)) then box_rad = 2.5
if not keyword_set(LINLIST) then $
   linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/lowd_ir_ohlines.lst' 

;; Determine slit width from header
slit_hdr   =  strtrim(sxpar(hdr, 'SLITNAME'))
slit_str = strmid(slit_hdr,3,7)
slit_arcsec = double(slit_str[0])
slit = slit_arcsec/0.143 
pkwdth = slit
TOLER = slit/2.0D
FWHM = slit

badpix = WHERE(finite(image) NE 1, nbad)
IF nbad NE 0 THEN image[badpix] = 0.0D

dims = size(image, /dimens)
nx = dims[0]
ny = dims[1]

; generate left and right edge of slits
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge

trace = (left_edge + right_edge)/2.0D
arc1d = fltarr(ny)
; Median filtering is more robust against cosmics
FOR j = 0L, ny-1L DO BEGIN
   left  = floor(trace[j] - BOX_RAD)
   right = ceil(trace[j] + BOX_RAD)
   sub_arc = image[left:right, j]
   djs_iterstat, sub_arc, median = median, sigrej = 2.0
   arc1d[j] = median
ENDFOR

;savefile='nirspec_Kband_36.02.sav'
;restore,savefile
;incalib=calib
;x_identify,arc1d,calib $
;           ,incalib=incalib $
;           , xsize = 1200, ysize = 600 $
;           , linelist = linlist $
;           , PKSIG = pksig, PKWDTH = pkwdth, TOLER = TOLER $
;           , /THIN, /FWEIGHT
;;           , mfitstr = calib, mspec = archive_arc ,mshift=shift $
;; Open line list
x_arclist, linlist, lines
;; what is cross-dispers position angle
disp_ang=float(sxpar(hdr, 'DISPPOS'))
case strtrim(filter,2) of
   'NIRSPEC-5': filtnm = 'filter5'
   else: filtnm = 'Kband'
endcase
prefix=getenv('LONGSLIT_DIR') + '/calib/linelists/nirspec_'+filtnm+'_*.sav'
archive_files=fileandpath(file_search(prefix,COUNT=COUNT))
IF COUNT EQ 0 THEN message,'Cannot find wavelngth archive'
archive_ang=fltarr(count)
FOR kk=0L,COUNT-1L DO BEGIN
   split=strsplit(archive_files[kk],'nirspec_'+filtnm+'_',/extract)
   archive_ang[kk]=float(split[0])
ENDFOR
diff=min(abs(archive_ang-disp_ang),imin)
print, 'nirspec_waveig: Reading  ',  archive_files[imin]
restore, getenv('LONGSLIT_DIR') + '/calib/linelists/' + archive_files[imin]
title_string=repstr(archive_files[imin],'.sav','')
title_string=repstr(title_string,'_','-')
;;  Zero out
lines.flg_plt = 0

;; Cross-correlate to get shift from archived arc
step = lindgen(2*MXSHFT) - MXSHFT 
pads = dblarr(2*MXSHFT)
corr = c_correlate([pads, (arc1d <  1.0e6), pads] $
                   , [pads, (archive_arc <  1.0e6), pads], step, /double)
xcen = long_find_nminima(-corr, nfind = 1, width = 5 $
                         , minsep = 1, ypeak = ypeak1, npeak = npeak)
xcen = double(xcen)      
shft = -step[long(round(xcen))]
;shft = -157L ;; JXP to get NIRSPEC-5 with 37.1 (01Aug2014)
print, shft, FORMAT = '(%"Wave Shift: %d pixels")' 

;   Re-identify lines
x_templarc, arc1d, lines, calib, MSK = msk $
            , SHFT = shft, ALL_PK = all_pk, PKWDTH = pkwdth, TOLER = TOLER $
            , FORDR = 9, PKSIG = pksig, FLG = flg_templ $
            , /THIN, /FWEIGHT
if flg_templ EQ -1 then begin
    print, 'x_fitarc: Insufficient lines for AUTO!!'
    stop
endif

;; Check the number of good lines
gdfit = where(lines.flg_plt EQ 1, ngd)
    
;; FIRST (crude) FIT
tmp_fit = {fitstrct}
copy_struct, calib, tmp_fit $
             , EXCEPT_TAGS = ['FFIT', 'NRM', 'RMS']
tmp_fit.flg_rej = 1 
tmp_fit.niter = 3 
tmp_fit.maxrej = 10 
tmp_fit.minpt = 5
tmp_fit.hsig = sigrej
tmp_fit.lsig = sigrej
tmp_fit.nord = (3 < calib.nord) ; Set first fit to 3
fin_fit = tmp_fit
fin_fit.nord = calib.nord

fit = x_fitrej(lines[gdfit].pix, lines[gdfit].wave, FITSTR = tmp_fit $
               , REJPT = rejpt)
if fit[0] EQ -1 then begin
    print, 'x_fitarc: AUTO Failed!!'
    stop
endif
    
;; Grab new lines
lines.flg_plt = 0
shft = 0L                       ; No shift expected now!
;; AutoID lines
x_templarc, arc1d, lines, tmp_fit, MSK = msk, PKWDTH = pkwdth, TOLER = TOLER $
            , FORDR = 9, SHFT = shft, PKSIG = pksig, /THIN, /FWEIGHT

;; Check the number of good lines
gdfit = where(lines.flg_plt EQ 1, ngd)
fin_fit = calib

fit = x_fitrej(lines[gdfit].pix, lines[gdfit].wave, FITSTR = fin_fit $
               , REJPT = rejpt, GDPT = gdpt)
ngd =  n_elements(gdpt)

; Output some QA
if keyword_set(QAFILE) then begin
    wv = x_calcfit(dindgen(ny), fitstr = fin_fit)
    dwv =  djs_median(abs(wv - shift(wv, 1)))
    rms = fin_fit.rms/dwv
    niri_waveqa, lines[gdfit], fit, wv, arc1d, rejpt, rms, QAFILE, title_string
    !p.multi = 0
ENDIF

pixset = long_wavepix(image, tset_slits, FWHM = FWHM, pkwdth = pkwdth $
                      , toler = toler, CHK = CHK)
piximg = long_wpix2image(pixset, tset_slits, XFIT = fin_fit, waveimg = waveimg)
; Apply heliocentric correction
ra = sxpar(hdr, 'RA')
dec = sxpar(hdr, 'DEC')
equinox = sxpar(hdr, 'EQUINOX')
jd = 2400000.5D + sxpar(hdr, 'MJD_OBS') 
helio = (-1.0d)*x_keckhelio(RA, DEC, EQUINOX, jd = JD)
hel_corr = sqrt( (1.d + helio/299792.458) / (1.d - helio/299792.458) )
print, 'nirspec_waveimg: Heliocentric correction applied -- ', helio, hel_corr
waveimg = waveimg*hel_corr

IF KEYWORD_SET(OUTFILE) THEN BEGIN
   ;;--------------                                
   ;;  Write output to wavefile
   splog, 'Writing output file'
   sxdelpar, hdr, 'NAXIS2'
   sxdelpar, hdr, 'NAXIS1'
   sxdelpar, hdr, 'NAXIS'
   sxaddpar, hdr, 'BITPIX', -64
   mwrfits, waveimg, outfile, hdr, /create 
   mwrfits, pixset, outfile 
   savefile = repstr(outfile, '.fits', '.sav')
   save,arc1d,fin_fit,file=savefile
ENDIF

RETURN, waveimg
END
