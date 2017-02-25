FUNCTION SOFI_WAVEIMG, image, tset_slits, hdr $
                       , piximg = piximg, MXSHIFT = MXSHIFT $
                       , SIGREJ = SIGREJ $
                       , BOX_RAD = BOX_RAD, QAFILE = QAFILE, CHK = CHK $
                       , CALIB = CALIB
  
;filename = '/Users/joe/DATA/SOFI_DATA/2011-09-22/SOFI_0220.fits'
;sofi_proc, filename, image, hdr = hdr
;qafile = 'test_H_sofi.ps'

if not keyword_set(MXSHFT) then MXSHFT = 15L
if not keyword_set(SIGREJ) then sigrej = 2.
if (NOT keyword_set(box_rad)) then box_rad = 5L
if not keyword_set(LINLIST) then $
   linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/oh_linelist.lst' 
;; Determine slit width from header
ind_slit = WHERE(stregex(hdr, 'HIERARCH ESO INS OPTI1 NAME', /bool))
slit_width = double(strmid(hdr[ind_slit], 42, 3))
plate_scale = 0.288d          
slit = slit_width/plate_scale ;; slit width in pixels
pkwdth = slit
toler = slit/2.0d > 2.0d
FWHM = slit

badpix = WHERE(finite(image) NE 1, nbad)
IF nbad NE 0 THEN image[badpix] = 0.0D

dims = size(image, /dimens)
nx = dims[0]
ny = dims[1]

; generate left and right edge of slits
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge

ind_setup = WHERE(stregex(hdr, 'HIERARCH ESO INS MODE', /bool))
setup = strcompress(repstr(strmid(hdr[ind_setup], 32, 13), "'", ' '), /rem)
CASE setup OF
   'LONG_SLIT_RED': BEGIN
      restore, file = getenv('XIDL_DIR') + $
               '/Spec/Longslit/calib/linelists/sofi_HK.sav'
      title_string = 'SOFI H+K grism'
      arc_offset = 0.0
   END
   'LONG_SLIT_H': BEGIN 
      restore, file = getenv('XIDL_DIR') + $
               '/Spec/Longslit/calib/linelists/sofi_H.sav'
      title_string = 'SOFI H grism'
      arc_offset = 30
   END
   'LONG_SLIT_K': BEGIN 
      splog, 'K is not yet supported!! WAVELENGTHS ARE JUNK!!!'
      pixset = long_wavepix(image, tset_slits, FWHM = FWHM, pkwdth = pkwdth $
                            , toler = toler, CHK = CHK)
      piximg = long_wpix2image(pixset, tset_slits)
      ;; Hack to create bogus wavelengths for now
      dlam = (2.52d4-2.0d4)/1024.0d
      waveimg = 2.0d4 + piximg*dlam
      RETURN, waveimg
      ;; 
      ;; Later use this code for the correct procedure
      ;restore, file = getenv('XIDL_DIR') + $
      ;         '/Spec/Longslit/calib/linelists/sofi_K.sav'
      ;title_string = 'SOFI K grism'
                                ;arc_offset = 0
   END
   ELSE: message, 'ERROR: Unknown setup'
ENDCASE

trace = (left_edge + right_edge)/2.0D + arc_offset
arc1d = fltarr(ny)
; Median filtering is more robust against cosmics
FOR j = 0L, ny-1L DO BEGIN
   left  = floor(trace[j] - BOX_RAD)
   right = ceil(trace[j] + BOX_RAD)
   sub_arc = image[left:right, j]
   djs_iterstat, sub_arc, median = median, sigrej = 2.0
   arc1d[j] = median
ENDFOR
;; Uncomment this if you want to dump an arc
;save, arc1d, file = 'sofi_arcH.sav'
;stop

;; Open line list
x_arclist, linlist, lines

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
print, shft, FORMAT = '(%"Wave Shift: %d pixels")' 

;   Re-identify lines
psig = 4.0D ;; this was default value
;; This is new code to improve the continuum fit of the arc which 
;; goes into the peak finding. Should be implemented everywhere that
;; uses x_templarc????
bkspace = 20.0
spec_set = bspline_iterfit(findgen(ny), arc1d $
                           , invvar = (arc1d GE 0.0 AND arc1d LE 1d6) $
                           , bkspace = bkspace $
                           , yfit = autofit $
                           , upper = upper, lower = lower, nord = 3 $
                           , maxrej = 10, outmask = outmask $
                           , /silent, /sticky)

x_templarc, arc1d, lines, calib, MSK = msk $
            , SHFT = shft, ALL_PK = all_pk, PKWDTH = pkwdth, TOLER = TOLER $
            , FORDR = FORDR, PKSIG = psig, FLG = flg_templ $
            , /THIN, /FWEIGHT, AUTOFIT = AUTOFIT
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
fsig = 4.0D 
x_templarc, arc1d, lines, tmp_fit, MSK = msk, PKWDTH = pkwdth, TOLER = TOLER $
            , FORDR = 9, SHFT = shft, PKSIG = fsig, /THIN, /FWEIGHT $
            , AUTOFIT = AUTOFIT

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
    niri_waveqa, lines[gdfit], fit, wv, arc1d, rejpt, rms, QAFILE $
                 , title_string, /log
    !p.multi = 0
ENDIF

;; !!!!!!FIX ME PLEASE!!!!!
;; ????? switch to tspec version but must fix a small bug first ????
pixset = long_wavepix(image, tset_slits, FWHM = FWHM, pkwdth = pkwdth $
                      , toler = toler, CHK = CHK)
piximg = long_wpix2image(pixset, tset_slits, XFIT = fin_fit, waveimg = waveimg)
; Apply heliocentric correction
ra = sxpar(hdr, 'RA')
dec = sxpar(hdr, 'DEC')
equinox = sxpar(hdr, 'EQUINOX')
jd = 2400000.5D + sxpar(hdr, 'MJD-OBS') 
helio = (-1.0d)*x_keckhelio(RA, DEC, EQUINOX, jd = JD)
hel_corr = sqrt( (1.d + helio/299792.458) / (1.d - helio/299792.458) )
print, 'sofi_waveimg: Heliocentric correction applied -- ', helio, hel_corr
waveimg = waveimg*hel_corr
waveimg = waveimg/1d4 ;; Wavelengths in micron

RETURN, waveimg
END
