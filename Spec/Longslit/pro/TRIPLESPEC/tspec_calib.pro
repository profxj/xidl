
;path = '/Users/joe/DATA/P200_DATA/2011418/'
;prefix = 'tspec110418_'
;darkind = 45L + lindgen(15)
;dark_str = string(darkind, FORMAT = '(I4.4)')

;filenames = path + prefix + dark_str + '.fits'
;superdark = niri_superdark(filenames, /TSPEC)

;; illumflatfiles
path = '/Users/joe/DATA/P200_DATA/2011416/'
prefix = 'tspec110416_'
flatind =  130L + lindgen(15)
flat_str = string(flatind, FORMAT = '(I4.4)')
illumflatfiles =  path + prefix + flat_str + '.fits'
;; pixflatfiles
path = '/Users/joe/DATA/P200_DATA/2011416/'
prefix = 'tspec110416_'
flatind =  12L + lindgen(10)
flat_str = string(flatind, FORMAT = '(I4.4)')
flatfiles =  path + prefix + flat_str + '.fits'
;; scifile
scifile = path + prefix + '0082.fits'


;; Find slits
orderfile = 'tspec-orders.fits'
tset_slits = tspec_findslits(flatfiles[0], orderfile, ksize = 7)


pixflatfile = 'tspec-pixflat.fits'
illumflatfile = 'tspec-illumflat.fits'
pixfile = 'tspec-piximg.fits'
;; Make Flat
tspec_makeflat, orderfile, flatfiles, illumflatfiles, scifile $
                , pixflatfile, illumflatfile, pixfile, /chk

END
