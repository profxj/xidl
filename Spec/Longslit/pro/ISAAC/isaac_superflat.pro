PRO isaac_superflat, flatfiles, darkfiles, objfile, pixflatfile $
                     , tempdir = tempdir, CHK = CHK
  

if (keyword_set(tempdir)) then $
   spawn, '\mkdir -p '+tempdir

nflat = n_elements(flatfiles)
ndark = n_elements(darkfiles)
IF nflat NE ndark THEN $
   message, 'Number of flat files and dark files do not agree'

tempflat = djs_filepath('tempflat-' + fileandpath(flatfiles) $
                        , root_dir = tempdir)
;; Process the FLAT-DARK sequence for Domeflats
FOR ii = 0L, nflat-1L DO BEGIN
   isaac_proc, flatfiles[ii], flat1, ivar_flat1, hdr = scihdr
   isaac_proc, darkfiles[ii], dark1, ivar_dark1, hdr = hdr
   ;; Average darks and flats
   sig2_flat = (ivar_flat1 GT 0)/(ivar_flat1 + (ivar_flat1 EQ 0))
   sig2_dark = (ivar_dark1 GT 0)/(ivar_dark1 + (ivar_dark1 EQ 0))
   sig2_diff = sig2_dark + sig2_flat
   ivar_diff = (ivar_dark1 GT 0)*(ivar_flat1 GT 0)*(sig2_diff GT 0.0)/(sig2_diff + (sig2_diff EQ 0))
   mwrfits, float(flat1-dark1), tempflat[ii], scihdr, /create
   mwrfits, float(ivar_diff), tempflat[ii]
ENDFOR

dims = size(flat1, /dim)
nx = dims[0]
ny = dims[1]
tset_slits = isaac_slitset(nx, ny)
slitmask = long_slits2mask(tset_slits)
slitfile = djs_filepath('slits-ISAAC.fits', root_dir = tempdir)
mwrfits, slitmask, slitfile, /create
mwrfits, tset_slits, slitfile

;; Create a wavefile for input to long_superflat
isaac_proc, objfile, arcimg, arcivar, hdr = hdr
scifile = djs_filepath('sci-' + fileandpath(scifile), root_dir = tempdir)
wavefile = repstr(scifile, 'sci', 'wave')
waveimg = isaac_waveimg(arcimg, arcivar, tset_slits, hdr, scifile $
                        , pixset = pixset)
mwrfits, waveimg, wavefile, /create 
mwrfits, pixset, wavefile 

use_pixel = [lonarr(nflat) + 1L]
use_illum = [lonarr(nflat)]
long_superflat, tempflat, pixflatfile, illumflatfile $
                , slitfile = slitfile, wavefile = wavefile $
                , use_illum = use_illum, use_pixel = use_pixel $
                , tempdir = tempdir, slitsamp = 20.0, CHK = CHK, /SOFI 


RETURN
END
