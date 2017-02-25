slitfile = 'slits-N20080407S0006.fits'
rawpath = '/data/gladders/gemini/GN-2008A-Q-25/raw/'
domefile = rawpath + 'N20080407S0006.fits'
;wavefile = 'wave-N20080407S0194.fits'
;scifile  = 'N20080407S0007.fits.gz'

;SDSSJ1446
;domefile = '../raw/N20080702S0057.fits'
;scifile = '../raw/N20080702S0058.fits'
;arcfile = '../raw/N20080702S0136.fits'
;wavefile = 'wave-N20080702S0136.fits'
;slitfile = 'slits-N20080702S0057.fits'

slitmask = mrdfits(slitfile, 0)
tset_slits = mrdfits(slitfile, 1)
ystruct = mrdfits(slitfile, 2)

nx = tset_slits[0].dims[0]
nmy = tset_slits[0].dims[1]
idim = size(tset_slits[0].COEFF, /dim)
nslit = idim[1]

traceset2xy, tset_slits[0], y_tr, xl_tr
traceset2xy, tset_slits[1], y_tr, xr_tr

gmos_trnxy1to3, reform(xl_tr, nmy*nslit), reform(y_tr, nmy*nslit) $
                , x3 = xl, y3 = yrow
gmos_trnxy1to3, reform(xr_tr, nmy*nslit), reform(y_tr, nmy*nslit) $
                , x3 = xr, y3 = yrow
xl = reform(xl, nmy, nslit)
xr = reform(xr, nmy, nslit)
yrow = reform(yrow, nmy, nslit)

long_proc, domefile, flat
slitmask = mrdfits(slitfile)
xatv, flat, min = -20.0, max = 1.0d4, wv = slitmask
;long_proc, scifile, sciimg
;xatv, sciimg, min = -20.0, max = 1.0d4, wv = slitmask
;flat = mrdfits(wavefile)

xatvplot, xl, yrow, psym = 1, color = 1
xatvplot, xr, yrow, psym = 1, color = 2

;PLOTSYM, 0, 1.0, THICK = 1.0, /FILL
;FOR islit = 0L, nslit-1L DO BEGIN
   ;iplot = where(yrow[*, islit] GE ystruct[islit].YSLIT_MIN AND $
   ;              yrow[*, islit] LE ystruct[islit].YSLIT_MAX, nip)
;   xatvplot, reform(left[iplot, islit], nip), reform(yrow[iplot, islit], nip)$
;             , psym = 8, color = 2
;   xatvplot, right[iplot, islit], yrow[iplot, islit], psym = 8, color = 1
;ENDFOR


END
