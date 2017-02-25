
;wvmnx = [1.40d, 2.50d]
;dlam = 0.880/2048.0d
;resolution = 2000.0d
;flgd = 0
;linefile = 'LUCI_modelsky_OH_linelist_200H+K.lst'
;outfile = 'LUCI_modelsky_OH_linelist_200H+K.fits'
;fnslit = 2.0d ;; assume 0.50 slit for model
;pkwdth = 1.3*fnslit
;toler = fnslit/2.0d > 2.0d
;thin = 1
;fweight = 0
;nsig = 5.0d
;nearir_modelsky_linelist, resolution, linefile, outfile $
 ;                         , wvmnx = wvmnx, dlam = dlam, flgd = flgd $
 ;                         , pkwdth = pkwdth, toler = toler, thin = thin $
 ;                         , fweight = fweight, nsig = nsig

reduxpath = '/Users/joe/LUCI_redux_0415/'
rawpath = reduxpath + 'Raw/'

darkfiles = rawpath + 'luci.20150413.00' + $
            strcompress(string(lindgen(10) + 50), /rem) + '.fits'
flatfiles = rawpath + 'luci.20150413.00' + $
            strcompress(string(lindgen(10) + 60), /rem) + '.fits'

slitfile = reduxpath + 'slits-' + fileandpath(flatfiles[0])
slitthresh = 0.3
luci_slitmask, flatfiles[0], slitfile, darkfile = darkfiles[0] $
               , minslit = minslit $
               , GMOSLONG = gmoslong $
               , peakthresh = slitthresh $
               , y1 = slity1, y2 = slity2 $
               , ksize = ksize, nfind = nfind $
               , verbose = verbose

slitmask = xmrdfits(slitfile, 0)
tset_slits = xmrdfits(slitfile, 1)

objfile = rawpath + 'luci.20150413.0078.fits'
luci_proc, objfile, wv_img, wv_ivar, hdr = hdr_wv
scifile = 'sci-' + fileandpath(objfile)
;wvchk = 1
;calib = 1
waveimg = luci_waveimg(wv_img, wv_ivar, tset_slits, hdr_wv, scifile $
                       , piximg = piximg, CHK = WVCHK, CALIB = CALIB)

luci_superflat, flatfiles, darkfiles, slitfile, objfile, pixflatfile $
                , tempdir = tempdir



END
