;domefile = 'N20080407S0018.fits.gz'
;scifile = 'N20080407S0016.fits.gz'
;arcfile = 'N20080407S0192.fits.gz'
;slitfile = 'slitmask.fits'

;; 0915+3826
domefile = 'N20080407S0006.fits.gz'
scifile = 'N20080407S0007.fits.gz'
arcfile = 'N20080407S0194.fits.gz'
wavefile = 'wave-N20080407S0194.fits'
slitfile = 'slits-N20080407S0006.fits'

;SDSSJ1446
;domefile = '../raw/N20080702S0057.fits'
;scifile = '../raw/N20080702S0058.fits'
;arcfile = '../raw/N20080702S0136.fits'
;wavefile = 'wave-N20080702S0136.fits'
;slitfile = 'slits-N20080702S0057.fits'

;; Create slitmask
gmos_slitmask, domefile, slitfile ;, /CHK
;; Create a waveimg
gmos_wavesolve, arcfile, slitfile, wavefile


END
