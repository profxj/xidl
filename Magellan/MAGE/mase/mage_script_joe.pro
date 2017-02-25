
;filenames = path+'mage0275.fits'

;path = '/Users/jhennawi/mage/nite1/'
;filenames = path + ['mage0076.fits', 'mage0078.fits']
;pixfile   = '/Users/jhennawi/mage/piximg-mage0085.fits'
;slitfile = '/Users/jhennawi/mage/mage_archive_orders.fits'
;wavefile = '/Users/jhennawi/mage/wavelength.fits'


;; trace the slits and generate a mask of the order numbers
;; This file is a sky flat for tracing the orders. Right now we are 
;; using the 2.0" slit because I don't know the skyflat for 0.85" slit
slitfile = '/Users/jhennawi/MAGE/orders-mage0192.fits'
traceflat = '/Users/jhennawi/MAGE/nite1/mage0192.fits'
IF NOT KEYWORD_SET(slitfile) THEN BEGIN
   mage_proc, traceflat, trcimg, hdr = hdr
   tset_slits = mage_traceorders(trcimg, chk = 1)
   slitmask = long_slits2mask(tset_slits)
   ordermask = 0*slitmask
   ordermask[WHERE(slitmask GT 0)] = -slitmask[WHERE(slitmask GT 0)] + 21L
   mwrfits, ordermask, slitfile, /create
   mwrfits, tset_slits, slitfile
ENDIF ELSE BEGIN
   ordermask = xmrdfits(slitfile, 0)
   tset_slits = xmrdfits(slitfile, 1)
   ;slitmask = long_slits2mask(tset_slits)
   ;ordermask = 0*slitmask
   ;ordermask[WHERE(slitmask GT 0)] = -slitmask[WHERE(slitmask GT 0)] + 21L
ENDELSE


IF NOT KEYWORD_SET(PIXFILE) THEN BEGIN $
   arcfile   = '/Users/jhennawi/MAGE/nite1/mage0085.fits'
   pixfile =  '/Users/jhennawi/MAGE/piximg-mage0085.fits'
   mage_proc, arcfile, arcimg, pixflatfile = pixflatfile
   pixset = long_wavepix(arcimg, tset_slits, fwhm = 3.0 $
                         , box_radius = 3.0 $
                         , sig_thresh = 3.0 $
                         , pkwdth = 5.0 $
                         , TOLER = 2.0 $
                         , CHK = 1)
   piximg = long_wpix2image(pixset, tset_slits)
   mwrfits, piximg, pixfile, /create
   mwrfits, pixset, pixfile
ENDIF ;;ELSE piximg = xmrdfits(pixfile, 0)

doflat = 0

if (doflat GT 0) then begin
   ;; This is the slit illumination flat created from a sky flat. We
   ;; are currently using the 2.0" slit until I find out which
   ;; files skyflats with 0.85" slit. We should probably just
   ;; set this to 1.0 for orders 20 and 6
   flatpath = '/Users/jhennawi/MAGE/nite1/'
   illumfiles =  flatpath + ['mage0192.fits']
   junkflat = mage_superflat(illumfiles, slitfile, pixfile $
                             , /skyflat, illum = illumflat)
   bluefiles  = flatpath + ['mage0194.fits', 'mage0195.fits']
   blueflat   = mage_superflat(bluefiles, slitfile, pixfile $
                               , order_vec = [20, 19, 18] $
                               , /skyflat)
   greenfiles = flatpath + ['mage0052.fits', 'mage0053.fits']
   greenflat  = mage_superflat(greenfiles, slitfile, pixfile $
                               , order_vec = [17, 16, 15, 14] $
                               , /skyflat)
   redfiles   = flatpath + ['mage0028.fits', 'mage0029.fits']
   redflat    = mage_superflat(redfiles, slitfile, pixfile $
                               , order_vec = [13, 12, 11, 10, 9, 8, 7, 6] $
                               , skyflat = 0)
   flat = blueflat * (ordermask LE 19 AND ordermask GE 18) + $
          greenflat * (ordermask LE 17 AND ordermask GE 14) + $
          redflat * (ordermask LE 13 AND ordermask GE 6)
   flat[where(ordermask EQ 20)] = 1.0
   xatv, flat, /block, min = 0.9, max = 1.1, /block
   xatv, illumflat, /block, min = 0.1, max = 2.0
endif
stop
pixflatfile = "mage_pixflat.fits"

;nfiles = n_elements(filenames)
;IF nfiles GT 0 THEN BEGIN
;   sciimg = 0
;   var_tot = 0
;   FOR ifile = 0L, nfiles-1L DO BEGIN
;      mage_proc, filenames[ifile], sciimg1, scivar1, pixflatfile=pixflatfile
;      sciimg = sciimg + sciimg1
;      var_tot = var_tot + 1.0D/(scivar1 + (scivar1 EQ 0))
;      var_tot = var_tot*(scivar1 GT 0)
;   ENDFOR
;   scivar = (var_tot GT 0)/(var_tot + (var_tot EQ 0))
;ENDIF ELSE  mage_proc, filenames, sciimg1, scivar, pixflatfile=pixflatfile
mage_proc, filenames, sciimg, scivar, pixflatfile = pixflatfile

ximg = long_slits2x(tset_slits, edgmask = edgmask)
waveimg = transpose(reverse(xmrdfits(wavefile, 0, /fscale), 1))
;waveimg = xmrdfits(wavefile, 0, /fscale)

bsp = 0.6
FWHM = 3.0
objstruct1 = long_objfind(sciimg, tset_slits = tset_slits, invvar = sciivar $
                          , skymask = skymask1, objmask = objmask1 $
                          , nperslit = 1L, peakthresh = reduxthresh $
                          , fwhm = FWHM, ISLIT = ISLIT)
skyimage = long_skysub(sciimg, sciivar, piximg, slitmask, skymask1, edgmask $
                       , bsp = bsp, ISLIT = ISLIT)
;objstruct = long_objfind(sciimg-skyimage, tset_slits = tset_slits $
;                         , invvar = sciivar, skymask = skymask $
;                         , objmask = objmask, nperslit = 1L $
;                         , peakthresh = reduxthresh $
;                         , fwhm = FWHM, ISLIT = ISLIT)


;splog, 'Redoing global sky subtraction'
;skyimage = long_skysub(sciimg, sciivar, piximg, slitmask, skymask, edgmask $
;                       , bsp = bsp, ISLIT = ISLIT)
; xatvplot, objstruct.xpos, objstruct.ypos, psym = 3
box_rad = 5L


velpix = 22.0d

;; The order struct is output by running mage_traceorders.pro
ordr_strct = xmrdfits('/home/simcoe/MagE/apr08/RS/redux1/OStr_new_mage.fits', 1)
obj_strct  = m_mkobjstr(15)

qafil = 'QAFntobj.ps'

;sciivar = 1.0d/(skyimage + (skyimage GT 0)) * (skyimage GT 0)

sciivar = 0.0*sciivar
sciivar[where(skyimage GT 0)] = 1.0d/(skyimage[where(skyimage GT 0)])

m_fntobj, obj_strct, ordr_strct, sciimg-skyimage, sciivar, qafil

xatv, (sciimg-skyimage)*(slitmask GT 0), min = -20.0, max = 200.0 $
      , wv = waveimg
for qq=0L,14 do xatvplot, obj_strct[qq].trace[0:2048-1], findgen(2048)

m_extechopt, sciimg, sciimg-skyimage, sciivar, ordr_strct, obj_strct, velpix, /chk, img_arc=alog10(waveimg), skyfil=skyimage, helio=0.0, obj_name="junk", ordermask=ordermask

;flux_box = extract_boxcar(sciimg-skyimage, objstruct.xpos, objstruct.ypos $
;                      , radius = box_rad)
;box_denom = extract_boxcar((waveimg GT 0.0), objstruct.xpos, objstruct.ypos $
;                           , radius = box_rad)
;wave_box = extract_boxcar(waveimg, objstruct.xpos, objstruct.ypos $
;                          , radius = box_rad)/(box_denom + (box_denom EQ 0))
;varimg = 1.0D/(scivar + (scivar EQ 0))
;var_box = extract_boxcar(varimg, objstruct.xpos, objstruct.ypos $
;                         , radius = box_rad)
;sky_box = extract_boxcar(skyimage, objstruct.xpos, objstruct.ypos $
;                         , radius = box_rad)
ord = 3
x_specplot, obj_strct[ord].fx, obj_strct[ord].sig, wav = obj_strct[ord].wave, inflg = 4



END
