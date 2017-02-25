;PRO nirspec_superflat, infiles, darkfiles $
;                       ,superdarkfile, pixflatfile, illumflatfile $
;                       , slitfile = slitfile $
;                       , objfile = objfile $
;                       , verbose = verbose, indir = indir $
;                       , tempdir = tempdir $
;                       , use_illum = use_illum $
;                       , use_pixel = use_pixel $
;                       , npoly = npoly, CHK = CHK $
;                       , _EXTRA = extra 


tempdir = 'Temp/'
;; LONG_SLIT_K, redo tonight
;;idom = [33, 39, 45]
;;itwi = [82] ;; image 84 is bad
;; LONG_SLIT_H
;idom = [3, 9, 15]
;itwi = [299, 301, 303]
;; LONG_SLIT_RED
path_dom = '/Users/joe/DATA/SOFI_DATA/2011-09-20/'
path_twi = '/Users/joe/DATA/SOFI_DATA/2011-09-22/'
;; need a science image for the wavelength solution
path_sci = '/Users/joe/DATA/SOFI_DATA/2011-09-22/'


idom = [4, 10, 16]
itwi = [85, 87]
;; Domeflats
FOR ii = 0L, n_elements(idom)-1L DO $
   IF ii EQ 0 THEN inds = [idom[ii] + lindgen(4)] ELSE inds = [inds, idom[ii] + lindgen(4)]
domefiles = path_dom + 'SOFI_' + string(inds, FORMAT = '(I4.4)') + '.fits'
ndome = n_elements(domefiles)/4L
tempdome = tempdir + 'TempDome-' + fileandpath(domefiles[4*lindgen(ndome)])
;; Twiflats
FOR ii = 0L, n_elements(itwi)-1L DO $
   IF ii EQ 0 THEN inds = [itwi[ii] + lindgen(2)] ELSE inds = [inds, itwi[ii] + lindgen(2)]
skyfiles = path_twi + 'SOFI_' + string(inds, FORMAT = '(I4.4)') + '.fits'
nsky = n_elements(skyfiles)/2L
temptwi = tempdir + 'TempTwi-' + fileandpath(skyfiles[2*lindgen(nsky)])

;; Create the slitmask file with one dark file and one flat file
darkfile = domefiles[0]
flatfile = domefiles[1]
slitfile = 'slitmask.fits'
luci_slitmask, flatfile, outfile, darkfile = darkfile, /SOFI
sofi_proc, flatfile, imag, darkfile = darkfile
slitmask = mrdfits(slitfile, 0)
tset_slits = mrdfits(slitfile, 1)

;; Create a superdark from the darks taken with the flats
fdarkfiles = strarr(ndome*2)
flatfiles = strarr(ndome*2)
ind = 0
FOR ii = 0L, ndome-1L DO BEGIN
   fdarkfiles[ind]   = domefiles[4*ii+0]
   flatfiles[ind]   = domefiles[4*ii+1]
   flatfiles[ind+1] = domefiles[4*ii+2]
   fdarkfiles[ind+1] = domefiles[4*ii+3]
   ind += 2
ENDFOR
flatdarkfile = 'flatdarkfile.fits'
darkimg = niri_superdark(fdarkfiles, outfile = fdarkfile, /SOFI)

;; identify the twilight flats and darks
tdarkfiles = strarr(nsky)
twifiles = strarr(nsky)
ind = 0
FOR ii = 0L, nsky-1L DO BEGIN
   twifiles[ind]     = skyfiles[2*ii+0]
   tdarkfiles[ind]   = skyfiles[2*ii+1]
   ind += 1
ENDFOR
twidarkfile = 'twiflatdarkfile.fits'
darkimg = niri_superdark(tdarkfiles, outfile = twidarkfile, /SOFI)

;; Create a dark file for cleaning up the wavelength solutions
darkfile = 'superdark120.fits'
darkfiles = path_sci + 'SOFI_' + string(lindgen(20) + 61, FORMAT = '(I4.4)') + '.fits'
darkimg = niri_superdark(darkfiles, outfile = darkfile, /SOFI)

;; Create a wavelength solution
scifile = path_sci + 'SOFI_0095.fits'
sofi_proc, scifile, arcimg, arcivar, hdr = hdr_arc, darkfile = darkfile
waveimg = luci_waveimg(arcimg, arcivar, tset_slits, hdr_arc $
                               , piximg = piximg, CHK = WVCHK, QAFILE = QAFILE)

nflat = n_elements(flatfiles)
ntwi = n_elements(twifiles)
infiles = [flatfiles, twifiles]
use_pixel = [lonarr(nflat) + 1L, lonarr(ntwi)]
use_illum = [lonarr(nflat), lonarr(ntwi) + 1L]
darkfiles = [replicate(flatdarkfile, nflat), replicate(twidarkfile, ntwi)]
stop


pixflatfile = 'pixflat-' + fileandpath(flatfiles[0])
illumflatfile = 'illumflat' + fileandpath(twifiles[0])


chk = 1
long_superflat, infiles, pixflatfile, illumflatfile $
                , waveimg = waveimg, piximg = piximg  $
                , slitfile = slitfile $
                , darkfiles = darkfiles $
                , use_illum = use_illum, use_pixel = use_pixel $
                , tempdir = tempdir, slitsamp = 5.0, CHK = CHK, /SOFI

sofi

END
