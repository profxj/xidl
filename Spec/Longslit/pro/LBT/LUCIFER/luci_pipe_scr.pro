
;11:33 Start calibration MPIA_Hennawi_revised
;L2-COSM_1001+0258: flat, arc #288-301
;L5-Lock_1033+5834: flat, arc #302-315


;;10:23 MPIA-A Hennawi L2-COSM_1001+0258
;;sky, src, slitview, verification, #122-125

;;10:38 Spectra started, #126-149, sky clear, seeing 0.8-1.4"

;;11:41 Telluric, acq, slitview, verification (#150-152), spectra #153-156
;;Offsets done manually

reduxpath = '/Users/joe/LUCI_redux_0415/'
path = reduxpath + 'Raw/'
darkfile = path + 'luci.20141204.0292.fits'
flatfile = path + 'luci.20141204.0293.fits'
;darkfile = 'Henn_20110103_0091.fits'
;flatfile = 'Henn_20110103_0092.fits'
slitfile = reduxpath + 'slits-' + fileandpath(flatfile)
slitthresh = 0.3
luci_slitmask, flatfile, slitfile, minslit = minslit $
               , GMOSLONG = gmoslong $
               , peakthresh = slitthresh $
               , y1 = slity1, y2 = slity2 $
               , ksize = ksize, nfind = nfind $
               , biasfile = darkfile, verbose = verbose






END
