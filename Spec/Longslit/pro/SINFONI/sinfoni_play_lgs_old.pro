
IF KEYWORD_SET(MAKE_MODEL) THEN BEGIN
   wvmnx = [1.9d, 2.5d] ;; expand the ragne to get slightly more coverage
   dlam = (2.45d - 1.95d)/2048.0d
   lam_cen = 1.95d + (2.45d - 1.95d)/2.0d
;; According to manual, assume resolution is actually 2 pixels
   resolution = lam_cen/(2.0*dlam)
   flgd = 0
   linefile = 'SINFONI_modelsky_OH_linelist_0.025_K.lst'
   outfile = 'SINFONI_modelsky_OH_linelist_0.025_K.fits'
   fnslit = 2.0d ;; slits are 2x plate scale
   pkwdth = 1.3*fnslit
   toler = fnslit/2.0d > 2.0d
   thin = 1
   fweight = 0
   nsig = 5.0d
   nearir_modelsky_linelist, resolution, linefile, outfile $
                             , wvmnx = wvmnx, dlam = dlam, flgd = flgd $
                             , pkwdth = pkwdth, toler = toler, thin = thin $
                             , fweight = fweight, nsig = nsig
ENDIF

path = '/Users/joe/SINFONI_K/data_with_raw_calibs/'
reduxpath = '/Users/joe/SINFONI_K/redux/'
flatfiles = ['SINFO.2014-07-25T12_35_42.142.fits' $
             , 'SINFO.2014-07-25T12_36_08.237.fits' $
             , 'SINFO.2014-07-25T12_39_02.749.fits']
flatdarkfiles = ['SINFO.2014-07-25T12_37_20.895.fits' $
                 , 'SINFO.2014-07-25T12_38_30.905.fits']
skyfiles = ['SINFO.2014-07-25T02_11_57.810.fits' $
            , 'SINFO.2014-07-25T03_45_53.093.fits']
objfiles = ['SINFO.2014-07-25T02_16_35.614.fits' $
            , 'SINFO.2014-07-25T03_36_51.556.fits']

darkfile = 'SINFO.2014-07-25T11_23_19.658.fits'

slitfile = 'slits-' + flatfiles[0]
superdarkfile = 'dark-' + darkfile
superflatdarkfile = 'flatdark-' + flatdarkfiles[0]
wavefile = 'wave-' + skyfiles[0]
savefile = repstr(wavefile, '.fits', '.sav')

IF file_test(superdarkfile) EQ 0 THEN BEGIN 
   darkimg = niri_superdark(darkfile, /SINFONI)
   mwrfits, darkimg, superdarkfile, /create
ENDIF

IF file_test(slitfile) EQ 0 THEN BEGIN 
   sinfoni_slitmask, path + flatfiles[0], slitfile $
                     , darkfile = superdarkfile  $
                     , minslit = minslit $
                     , peakthresh = slitthresh $
                     , y1 = slity1, y2 = slity2 $
                     , ksize = ksize, nfind = nfind $
                     , verbose = verbose
   slits = mrdfits(slitfile, 0)
   sinfoni_proc, flatfile, flatimg
;xatv, flatimg, /align
;xatv, slits, /align, max = 1.0
ENDIF


;sinfoni_proc, skyfile, arcimg, arcivar $
;              , hdr = hdr_arc, darkfile = superdarkfile, verbose = verbose

tset_slits = mrdfits(slitfile, 1)
IF file_test(wavefile) EQ 0 THEN BEGIN
   ;;wvchk = 1
   QAFILE = 'SINFONI_wave_QA.ps'
   waveimg = sinfoni_waveimg(skyfiles, tset_slits, hdr_arc $
                             , piximg = piximg, CHK = WVCHK, QAFILE = QAFILE $
                             , darkfile = superdarkfile $
                             , outfile = wavefile, savefile = savefile)
ENDIF

nflat = n_elements(flatfiles)
use_pixel = [lonarr(nflat) + 1L]
use_illum = [lonarr(nflat) + 1L]
pixflatfile = 'pixflat-' + flatfiles[0]
illumflatfile = 'illumflat-' + flatfiles[0]
;;chk = 1
;; Illumination function fits in long_flatfield_specillum.pro
;; appear to be choking for a few slits and returning zero. This needs
;; to be debugged.
IF file_test(pixflatfile) EQ 0 THEN BEGIN
   IF file_test(superflatdarkfile) EQ 0 THEN BEGIN 
      darkimg = niri_superdark(flatdarkfiles, /SINFONI)
      mwrfits, darkimg, superflatdarkfile, /create
   ENDIF
   tempdir = 'Temp/'
   long_superflat, flatfiles, pixflatfile, illumflatfile $
                   , slitfile = slitfile, wavefile = wavefile $
                   , darkfile = superflatdarkfile $
                   , use_illum = use_illum, use_pixel = use_pixel $
                   , tempdir = tempdir, slitsamp = 5.0, CHK = CHK, /SINFONI
;;sinfoni_proc, flatfiles[0], flat, illum = illumflatfile, pixflat = pixflatfile
;;sinfoni_proc, flatfiles[0], flat0
;xatv, flat, min = -20.0, max = 40000, /align
;xatv, flat0, min = -20.0, max = 40000, /align
ENDIF
stop
;illumflatfile = 'Temp/tempillum-SINFO.2014-07-25T12_35_42.142.fits'
;pixflatfile = 'Temp/temppixel-SINFO.2014-07-25T12_35_42.142.fits'

sinfoni_proc, objfiles[0], obj_img, ivar_obj, hdr = hdr_obj $
              , pixflat = pixflatfile, illumflat = illumflatfile
sinfoni_proc, skyfiles[0], sky_img, ivar_sky, hdr = hdr_sky $
              , pixflat = pixflatfile, illumflat = illumflatfile
;; Do first sky-subtraction
mask_obj = (ivar_obj GT 0.0)
mask_sky = (ivar_sky GT 0.0)
var_obj = (ivar_obj GT 0.0)/(ivar_obj + (ivar_obj LE 0.0))
var_sky = (ivar_sky GT 0.0)/(ivar_sky + (ivar_sky LE 0.0))
obj_min_sky = obj_img - sky_img
var_obj_min_sky = var_obj + var_sky
mask = (mask_obj GT 0.0) AND (mask_sky GT 0.0)
ivar = float(mask)/(var_obj_min_sky + (var_obj_min_sky EQ 0.0))



END
