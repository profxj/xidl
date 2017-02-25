filename = '/Users/joe/luci_data/luci.20100312.0074.fits'
obj = mrdfits(filename, 0, hdr, /fscale)
sciimg = reverse(transpose(obj), 2)
gain = float(sxpar(hdr, 'ELECGAIN'))
rdnoise = float(sxpar(hdr, 'ENOISE'))
sciimg = gain*sciimg
sciivar = 1.0/(abs(sciimg - sqrt(2.0)*rdnoise) +rdnoise^2)

slitfile = '/Users/joe/luci_data/slits-luci.20100312.0098.fits'
tset_slits = mrdfits(slitfile, 1)
xshift = long_xcorr_slits(sciimg, tset_slits, /shift)
traceset2xy, tset_slits[0], rows, left_edge
traceset2xy, tset_slits[1], rows, right_edge
slitmask = long_slits2mask(tset_slits)
ximg = long_slits2x(tset_slits, TOL_EDG = 7.0, EDGMASK = EDGMASK)

;; Create pixel wavelength map
slit_arcsec = 0.8D
slit = slit_arcsec/0.25D ;; plate scale. Slit is slit-width in pixesl
pkwdth = slit
TOLER = slit/2.0D
FWHM = slit

;stop
;CHK = 1
pixset = long_wavepix(sciimg, tset_slits, FWHM = FWHM, pkwdth = pkwdth $
                      , toler = toler, CHK = CHK)
piximg = long_wpix2image(pixset, tset_slits, XFIT = fin_fit, waveimg = waveimg)

bsp = 0.7D
skyimage = long_skysub(sciimg, sciivar, piximg, slitmask, skymask1, edgmask $
                        , bsp = bsp)


END
