;+
; NAME:
;   isaac_proc
;
; PURPOSE:
;   Process an ISAAC near-IR image
;
; CALLING SEQUENCE:
;   isaac_proc, filename, [ flux, invvar, hdr= $
;    , flatfile=, adderr=, /verbose ]
;
; INPUTS:
;   filename   -  name of a image file
;
; OPTIONAL INPUTS:
;   flatfile- File with pixel flat
;   adderr     - Additional error to add to the formal errors, as a fraction
;                of the flux; default to 0.01 (1 per cent).
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   flux       - Image 
;   invvar     - Inverse variance image
;   hdr        - FITS header
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES CALLED:
;   divideflat
;   headfits()
;   mrdfits()
;   niri_oscan
;
; REVISION HISTORY:
;   6-June-2006  Written by J. Hennawi (UCB)
;-
;------------------------------------------------------------------------------

pro isaac_proc, filename, imag, invvar, hdr = hdr $
                , pixflatfile = pixflatfile, darkfile = darkfile $
                , illumflatfile = illumflatfile $
                , adderr = adderr, verbose = verbose, bin = bin $
                , gain = gain, rnoise = rdnoise, image = imag1 $
                , OSCAN_SUB = OSCAN_SUB
 
minval = 0.5
if (n_elements(adderr) EQ 0) then adderr = 0.01

hdr = xheadfits(filename)

if (size(hdr, /tname) NE 'STRING') then begin
    splog, 'Invalid FITS header for file ', filename
    imag = 0
    invvar = 0
    return
endif

;;  If FLUX and INVVAR are not requested, then return with just the header
if (arg_present(imag) EQ 0 AND arg_present(invvar) EQ 0) then return

;; ISAAC Hawaii Hg:Cd:Te 1024x1024 detector
GAIN = 4.5 ;; e-/ADU
RDNOISE = 11.0d ;; measure this

filelist = lookforgzip(filename)
if filelist[0] EQ '' then begin
    print, 'Could not find file named ', filename
    return
endif
hdr0 = xheadfits(filelist[0])
; lower wavelengths are up by default, so reverse to flip
imag1 = transpose(xmrdfits(filelist[0], 0, /silent))
;;imag = reverse(xmrdfits(filelist[0], 0, /silent), 2)


if size(hdr0, /tname) EQ 'INT' then begin
    print, 'Having trouble with ', filename
    return
endif

dims = size(imag1, /dim)
nx = dims[0]
ny = dims[1]

IF KEYWORD_SET(OSCAN_SUB) OR arg_present(invvar) THEN BEGIN
;; Compute the oscan which is needed for the 
;; This is a kludge for now. Replace this with a superdark to get
;; the noise statistics right. 
   l_mask = lonarr(nx, ny)
   l_mask[30:60, 65:420] = 1L
   l_mask[30:60, 620:900] = 1L
   r_mask = lonarr(nx, ny)
   r_mask[960:1000, 10:400] = 1L
   r_mask[960:1000, 600:830] = 1L
   il = WHERE(l_mask)
   ir = WHERE(r_mask)
   djs_iterstat, imag1[il], median = oscan_l, mean = mean_l, sigma = sig_l
   djs_iterstat, imag1[ir], median = oscan_r, mean = mean_r, sigma = sig_r
   oscan = fltarr(nx, ny)
   oscan[0:511, *] = oscan_l
   oscan[512:*, *] = oscan_r
ENDIF

;; Note we difference science frames, and subtract simultaneous darks
;; from the flats, so normally there is no need to do subtract the
;; 'oscan', which is actually a complex spatially varying dark current
;; amplifier/readout pattern. 
IF KEYWORD_SET(OSCAN_SUB) THEN imag = GAIN*(imag1-oscan) $
ELSE imag = GAIN*imag1

;; Always subtract off the oscan/dark for the inverse variance, since
;; otherwise we get the photon counting statistics wrong. Does the
;; dark current/pattern vary with time?? This is oscan is a kludge for
;; now. 
IF arg_present(invvar) THEN BEGIN
   imag_min_oscan = GAIN*(imag1-oscan)
   invvar = 1.0/(abs(imag_min_oscan - sqrt(2.0)*rdnoise) +rdnoise^2)
ENDIF
;; Add bad pixel mask here
if (arg_present(invvar) AND adderr GT 0) then begin
   gmask = invvar GT 0          ; =1 for good points
   invvar = gmask / (1.0 / (invvar + (1-gmask)) + $
                     adderr^2 * (abs(imag_min_oscan))^2)
endif

ndim = size(imag, /n_dimen)
dims = size(imag, /dimens)

if (keyword_set(darkfile)) then begin
   darkimg = xmrdfits(darkfile, 0, biashdr, $
                      silent = (keyword_set(verbose) EQ 0))
   if (size(darkimg, /n_dimen) NE ndim $
       OR total(size(darkimg, /dimens) NE dims) NE 0) then $
          message, 'Dimensions of image and darkimg image do not agree'
   imag = imag - darkimg
endif

if (keyword_set(pixflatfile)) then begin
   flatimg = xmrdfits(pixflatfile, 0, biashdr, $
                      silent = (keyword_set(verbose) EQ 0))
   if (size(flatimg, /n_dimen) NE ndim $
       OR total(size(flatimg, /dimens) NE dims) NE 0) then $
          message, 'Dimensions of image and flat image do not agree'
;   Don't mask the flat
   divideflat, imag, flatimg, invvar = invvar, minval = 0.0 $
               , quiet = (keyword_set(verbose) EQ 0)
;   This will act like a bad pixel mask
   IF ARG_PRESENT(invvar) THEN BEGIN
      badmin = 0.7
      badmax = 1.3
      badpixmask = (flatimg GT badmin AND flatimg LT badmax)
      invvar = invvar*badpixmask
   ENDIF
ENDIF

if (keyword_set(illumflatfile)) then begin
   illumflat = xmrdfits(illumflatfile, 0, biashdr, $
                        silent = (keyword_set(verbose) EQ 0))
   if (size(illumflat, /n_dimen) NE ndim $
       OR total(size(illumflat, /dimens) NE dims) NE 0) then $
          message, 'Dimensions of image and flat image do not agree'
   ;; Don't apply illumination corrections larger than 30%
   gdpix = WHERE(illumflat GT 0.6D AND illumflat LT 1.4)
   imag[gdpix]  = imag[gdpix]/illumflat[gdpix]
   IF ARG_PRESENT(invvar) THEN $
      invvar[gdpix] = invvar[gdpix]*illumflat[gdpix]^2 
ENDIF

return
end

;------------------------------------------------------------------------------
