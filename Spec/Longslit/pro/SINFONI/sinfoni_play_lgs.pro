
;; This generates the model spectrum of IDing the sky-lines
IF KEYWORD_SET(MAKE_MODEL) THEN BEGIN
   wvmnx = [1.9d, 2.5d] ;; expand the ragne to get slightly more coverage
   dlam = (2.45d - 1.95d)/2048.0d
   lam_cen = 1.95d + (2.45d - 1.95d)/2.0d
;; According to manual, assume resolution is actually 2 pixels
   resolution = lam_cen/(2.0*dlam)
   flgd = 0
   linefile = 'SINFONI_modelsky_OH_linelist_0.25_K.lst'
   outfile = 'SINFONI_modelsky_OH_linelist_0.25_K.fits'
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

path = '/Users/joe/SINFONI_082015/2015-08-04/'
reduxpath = '/Users/joe/SINFONI_082015/redux_lgs/'
flatfiles = path + ['SINFO.2015-08-02T11_19_13.232.fits'   $
                    , 'SINFO.2015-08-02T11_19_38.440.fits' $ 
                    , 'SINFO.2015-08-02T11_22_34.753.fits' $     
                    , 'SINFO.2015-08-02T11_22_59.078.fits' $
                    , 'SINFO.2015-08-02T11_25_54.478.fits']
flatdarkfiles = path + ['SINFO.2015-08-02T11_18_44.418.fits' $
                        , 'SINFO.2015-08-02T11_20_53.819.fits' $
                        , 'SINFO.2015-08-02T11_22_02.912.fits' $
                        , 'SINFO.2015-08-02T11_24_10.856.fits' $
                        , 'SINFO.2015-08-02T11_25_22.640.fits']

skyfiles = path + 'SINFONI_IFS_SKY217_' + $
           string(lindgen(4) + 1L, format = '(I4.4)') + '.fits'

objfiles = path + 'SINFONI_IFS_OBS217_' + $
           string(lindgen(8) + 1L, format = '(I4.4)') + '.fits'

darkfiles = path + [ 'SINFO.2015-08-02T10_50_13.961.fits' $
                     , 'SINFO.2015-08-02T10_55_34.430.fits' $
                     , 'SINFO.2015-08-03T10_51_17.548.fits' $ 
                     , 'SINFO.2015-08-03T10_56_39.674.fits' $ 
                     , 'SINFO.2015-08-03T11_01_57.673.fits' $  
                     , 'SINFO.2015-08-04T11_09_53.539.fits' $ 
                     , 'SINFO.2015-08-04T11_15_15.672.fits' $   
                     , 'SINFO.2015-08-04T11_20_34.486.fits'] 

slitfile = reduxpath + 'slits-' + fileandpath(flatfiles[0])
superdarkfile = reduxpath + 'dark-' + fileandpath(darkfiles[0])
superflatdarkfile = reduxpath + 'flatdark-' + fileandpath(flatdarkfiles[0])

wavefile = reduxpath + 'wave-' + fileandpath(skyfiles[0])
savefile = repstr(wavefile, '.fits', '.sav')

IF file_test(superflatdarkfile) EQ 0 THEN BEGIN 
   darkimg = niri_superdark(flatdarkfiles, /SINFONI)
   mwrfits, darkimg, superflatdarkfile, /create
ENDIF

stop
;IF file_test(slitfile) EQ 0 THEN BEGIN
;   slit_edg_twk = 1
   sinfoni_slitmask, flatfiles[0], flatdarkfiles, slitfile $
                     , minslit = minslit $
                     , peakthresh = slitthresh $
                     , y1 = slity1, y2 = slity2 $
                     , ksize = ksize, nfind = nfind $
                     , verbose = verbose ;;, slit_edg_twk = slit_edg_twk
   stop
   slitmask = mrdfits(slitfile, 0)
   sinfoni_proc, flatfiles[0], flatimg
;xatv, flatimg, /align
;xatv, slits, /align, max = 1.0
;ENDIF


IF file_test(superdarkfile) EQ 0 THEN BEGIN 
   darkimg = niri_superdark(darkfiles, /SINFONI)
   mwrfits, darkimg, superdarkfile, /create
ENDIF

slitmask = mrdfits(slitfile, 0)
tset_slits = mrdfits(slitfile, 1)
IF file_test(wavefile) EQ 0 THEN BEGIN
   ;;wvchk = 1
   QAFILE = repstr(wavefile, '.fits', '.ps')
   waveimg = sinfoni_waveimg(skyfiles, tset_slits $
                             , piximg = piximg, CHK = WVCHK, QAFILE = QAFILE $
                             , darkfile = superdarkfile $
                             , outfile = wavefile, savefile = savefile)
ENDIF



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
pixflatfile = reduxpath + 'pixflat-' + fileandpath(flatfiles[0])
illumflatfile = reduxpath + 'illumflat-' + fileandpath(flatfiles[0])
;;chk = 1
;; Illumination function fits in long_flatfield_specillum.pro
;; appear to be choking for a few slits and returning zero. This needs
;; to be debugged.
IF file_test(pixflatfile) EQ 0 THEN BEGIN
   tempdir = reduxpath + 'Temp/'
   ;;chk = 1
   long_superflat, fileandpath(flatfiles), fileandpath(pixflatfile), fileandpath(illumflatfile) $
                   , slitfile = slitfile, wavefile = wavefile $
                   , darkfile = superflatdarkfile $
                   , use_illum = use_illum, use_pixel = use_pixel $
                   , tempdir = tempdir, slitsamp = 5.0, CHK = CHK, /SINFONI
;;sinfoni_proc, flatfiles[0], flat, illum = illumflatfile, pixflat = pixflatfile,dark=superdarkfile
;;sinfoni_proc, flatfiles[0], flat0,dark=superdarkfile
;xatv, flat, min = -20.0, max = 40000, /align
;xatv, flat0, min = -20.0, max = 40000, /align
ENDIF

wave_qa_file = repstr(reduxpath + $
                      'wave-' + fileandpath(skyfiles[0]), '.fits', '.ps')
scifile = reduxpath + 'sci-' + fileandpath(objfiles[4])

;; Compute the wavelength image for this source from the skyfile stack
waveimg = sinfoni_waveimg(skyfiles, tset_slits $
                          , piximg = piximg, CHK = WVCHK, QAFILE = wave_qa_file $
                          , darkfile = darkfile $
                          , pixflatfile = pixflatfile $
                          , illumflatfile = illumflatfile $
                          , outfile = outfile, savefile = savefile)

chk = 1
obj_min_sky = sinfoni_skysub(objfiles[4], skyfiles[2], tset_slits, piximg $
                             , ivar = ivar $
                             , sky_resids = sky_resids $
                             , sky_model = sky_model $
                             , hdr_obj = scihdr $
                             , plate_scale = plate_scale $
                             , darkfile = superdarkfile $
                             , pixflatfile = pixflatfile $
                             , illumflatfile = illumflatfile, chk = chk)
IF KEYWORD_SET(FWHM1) THEN FWHM = FWHM1 ELSE FWHM = 4.5
objstruct = long_objfind(obj_min_sky - sky_resids, tset_slits = tset_slits $
                         , FWHM = FWHM, OBJMASK = OBJMASK_POS $
                         , NPERSLIT = KEYWORD_SET(TELLURIC) $
                         , OBJTHRESH = KEYWORD_SET(TELLURIC)*0.1 $
                         , peakthresh = peakthresh)
final_struct = 0
dimt = size(tset_slits.coeff, /dimen)
nslit = dimt[1]
nobj = n_elements(objstruct)
;----------
; Loop over each object and extract
FOR iobj = 0L, nobj-1L DO BEGIN
   thismask = (slitmask EQ objstruct[iobj].SLITID)
   extract_obj = gnirs_extract(obj_min_sky-sky_resids, ivar, waveimg $
                               , thismask, sky_model, objstruct[iobj] $
                               , plate_scale, TELLURIC = TELLURIC $
                               , SN_GAUSS = SN_GAUSS)
   final_struct = struct_append(final_struct, extract_obj)
ENDFOR

;----------
; Write output file

splog, 'Writing FITS file ', scifile
mwrfits, float(obj_min_sky), scifile, scihdr, /create
mwrfits, float(sky_resids), scifile
mwrfits, float(ivar), scifile
mwrfits, float(waveimg), scifile
mwrfits, final_struct, scifile



END
