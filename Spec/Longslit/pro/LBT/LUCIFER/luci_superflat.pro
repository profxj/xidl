;PRO luci_superflat, flatfiles, darkfiles, twifiles, twidark $
;                    , slitfile, objfile, pixflatfile $
;                    , tempdir = tempdir, CHK = CHK


reduxpath = '/Users/joe/LUCI_redux_0415/'
slitfile = reduxpath + 'slits-luci.20150413.0060.fits'

nflat = 5
;nflat = 2L
tempdir = 'Temp/'
path = '/Users/joe/LUCI_data_april15/20150414/'
darkfiles = path + 'luci.20150414.00' + $
            strcompress(string(lindgen(nflat) + 37), /rem) + '.fits'
flatfiles = path + 'luci.20150414.00' + $
            strcompress(string(lindgen(nflat) + 47), /rem) + '.fits'
path2 = '/Users/joe/LUCI_data_april15/20150413/'
twifiles = path2 + 'luci.20150413.0' + $
           strcompress(string(lindgen(nflat) + 149), /rem) + '.fits'
twidark = path2 + 'luci.20150413.0' + $
          strcompress(string(lindgen(nflat) + 128), /rem) + '.fits'

if (keyword_set(tempdir)) then $
   spawn, '\mkdir -p '+tempdir


objfile = path2 + 'luci.20150413.0078.fits'

tset_slits = xmrdfits(slitfile, 1)
slitmask = long_slits2mask(tset_slits)
;; Create a wavefile for input to long_superflat
luci_proc, objfile, arcimg, arcivar, hdr = hdr_wv
scifile = djs_filepath('sci-' + fileandpath(objfile), root_dir = tempdir)
wavefile = repstr(scifile, 'sci', 'wave')
waveimg = luci_waveimg(arcimg, arcivar, tset_slits, hdr_wv, scifile $
                       , pixset = pixset)
mwrfits, waveimg, wavefile, /create 
mwrfits, pixset, wavefile 


nflat = n_elements(flatfiles)
ndark = n_elements(darkfiles)
ntwi = n_elements(twifiles)
ntwidark = n_elements(twidark)
IF nflat NE ndark THEN $
   message, 'Number of flat files and dark files do not agree'

tempflat = djs_filepath('tempflat-' + fileandpath(flatfiles) $
                        , root_dir = tempdir)
pixflatfile = 'pixflat.fits'
illumflatfile = 'illumflat.fits'

itell = WHERE((waveimg GE 1.340D AND waveimg LE 1.5D) OR $
              (waveimg GE 1.785D AND waveimg LE 2.0D), ntell)
;; Process the FLAT-DARK sequence for Domeflats
FOR ii = 0L, nflat-1L DO BEGIN
   luci_proc, flatfiles[ii], flat1, ivar_flat1, hdr = scihdr
   luci_proc, darkfiles[ii], dark1, ivar_dark1, hdr = hdr
   ;; Average darks and flats
   sig2_flat = (ivar_flat1 GT 0)/(ivar_flat1 + (ivar_flat1 EQ 0))
   sig2_dark = (ivar_dark1 GT 0)/(ivar_dark1 + (ivar_dark1 EQ 0))
   sig2_diff = sig2_dark + sig2_flat
   ivar_diff = (ivar_dark1 GT 0)*(ivar_flat1 GT 0)*(sig2_diff GT 0.0)/(sig2_diff + (sig2_diff EQ 0))
   ;; Mask telluric absorption region from twiflat
   mwrfits, float(flat1-dark1), tempflat[ii], scihdr, /create
   mwrfits, float(ivar_diff), tempflat[ii]
ENDFOR

temptwi = djs_filepath('temptwi-' + fileandpath(twifiles) $
                       , root_dir = tempdir)

;; Process the TWIFLAT and TWIDARK sequence for twilight flats
FOR ii = 0L, nflat-1L DO BEGIN
   luci_proc, twifiles[ii], flat1, ivar_flat1, hdr = scihdr
   luci_proc, twidark[ii], dark1, ivar_dark1, hdr = hdr
   ;; Average darks and flats
   sig2_flat = (ivar_flat1 GT 0)/(ivar_flat1 + (ivar_flat1 EQ 0))
   sig2_dark = (ivar_dark1 GT 0)/(ivar_dark1 + (ivar_dark1 EQ 0))
   sig2_diff = sig2_dark + sig2_flat
   ivar_diff = (ivar_dark1 GT 0)*(ivar_flat1 GT 0)*(sig2_diff GT 0.0)/(sig2_diff + (sig2_diff EQ 0))
   IF ntell GT 0 THEN ivar_diff[itell] = 0
   mwrfits, float(flat1-dark1), temptwi[ii], scihdr, /create
   mwrfits, float(ivar_diff), temptwi[ii]
ENDFOR

tempfiles = [tempflat, temptwi]
use_pixel = [lonarr(nflat) + 1L, lonarr(ntwi)]
use_illum = [lonarr(nflat), lonarr(ntwi) + 1L]


long_superflat, tempfiles, pixflatfile, illumflatfile $
                , slitfile = slitfile, wavefile = wavefile $
                , use_illum = use_illum, use_pixel = use_pixel $
                , tempdir = tempdir, slitsamp = 5.0, CHK = CHK, /LUCI


RETURN
END
