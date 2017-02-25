;+
; NAME:
;   tspec_superdark
;
; PURPOSE:
;   Used to generate a Super dark file for TRIPLESPEC
;
; CALLING SEQUENCE:
;   tspec_superdark, ['files']
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   tset_slits - 2-element array of trace sets, where the first defines
;                the starting slit positions, and the second one defines
;                the ending slit positions
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   ;
; REVISION HISTORY:
;   
;-
;------------------------------------------------------------------------------
PRO TSPEC_SUPERDARK, darklist, superdarkfile=superdarkfile, badpixmaxkfile=badpixmaskfile, $
                     SLITFILE=slitfile, BAD_VALUE=bad_value


readcol, darklist, darkfiles, format='A'
if not keyword_set(SUPERDARKFILE) then SUPERDARKFILE = 'tspec-superdark.fits'
if not keyword_set(BADPIXMASKFILE) then BADPIXMASKFILE = 'tspec-badpix.fits'
if not keyword_set(BAD_VALUE) then BAD_VALUE = 10.

;; First make the superdark image
if x_chkfil(superdarkfile) then darkimg = xmrdfits(superdarkfile) else begin
   darkimg = niri_superdark(darkfiles, /TSPEC)
   mwrfits, darkimg, superdarkfile, /create
endelse

dims = size(darkimg, /dim)
nx = dims[0]
ny = dims[1]
ndark = n_elements(darkfiles)

;; Generate the slit
if not keyword_set(SLITFILE) then begin
   tset_slits = kast_slitset(nx, ny)
   slitmask = long_slits2mask(tset_slits)
   slitfile = 'slits-TSPEC_DARK.fits'
   mwrfits, slitmask, slitfile, /create
   mwrfits, tset_slits, slitfile
endif

;; output files
pixflatfile = 'pixdark-' + fileandpath(darkfiles[0])
illumflatfile = 'illumdark-' + fileandpath(darkfiles[0])

use_pixel = [lonarr(ndark) + 1L]
use_illum = [lonarr(ndark) + 1L]

if not x_chkfil(pixflatfile) and not keyword_set(REDO) then begin
   long_superflat, darkfiles, pixflatfile, illumflatfile $
                   , slitfile = slitfile $
                   , use_illum = use_illum, use_pixel = use_pixel $
                   , tempdir = tempdir, slitsamp = 20.0, /TSPEC ;;, /CHK ;;, /CHK
endif

normdark = xmrdfits(pixflatfile, 0)
badpixmask = lonarr(nx, ny) + 1L
badpix = WHERE(normdark GE BAD_VALUE, nbad)
IF nbad GT 0 THEN badpixmask[badpix] = 0L

xatv, badpixmask, max=2., /bloc

;; Write
mwrfits, badpixmask, badpixmaskfile, /create





;; The normalized pixel flat field image produced from the dark has
;; the spectral and spatial illumination pattern fit out. Cut this for
;; outlying hot pixels and use that to define the bad pixel mask. 

END
