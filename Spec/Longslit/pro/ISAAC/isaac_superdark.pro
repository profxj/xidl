PRO ISAAC_SUPERDARK, darkfiles, superdarkfile, badpixmaskfile


;; First make the superdark image
darkimg = niri_superdark(darkfiles, /ISAAC)
mwrfits, darkimg, superdarkfile, /create

ndark = n_elements(darkfiles)

dims = size(darkimg, /dim)
nx = dims[0]
ny = dims[1]
tset_slits = isaac_slitset(nx, ny)
slitmask = long_slits2mask(tset_slits)
slitfile = 'slits-ISAAC.fits'
mwrfits, slitmask, slitfile, /create
mwrfits, tset_slits, slitfile

;; output files
pixflatfile = 'pixdark-' + fileandpath(darkfiles[0])
illumflatfile = 'illumdark-' + fileandpath(darkfiles[0])

use_pixel = [lonarr(ndark) + 1L]
use_illum = [lonarr(ndark) + 1L]
long_superflat, darkfiles, pixflatfile, illumflatfile $
                , slitfile = slitfile $
                , use_illum = use_illum, use_pixel = use_pixel $
                , tempdir = tempdir, slitsamp = 20.0, /ISAAC ;;, /CHK ;;, /CHK


normdark = xmrdfits(pixflatfile, 0)
badpixmask = lonarr(nx, ny) + 1L
badpix = WHERE(normdark GE 3.0, nbad)
IF nbad GT 0 THEN badpixmask[badpix] = 0L
mwrfits, badpixmask, badpixmaskfile, /create





;; The normalized pixel flat field image produced from the dark has
;; the spectral and spatial illumination pattern fit out. Cut this for
;; outlying hot pixels and use that to define the bad pixel mask. 

END
