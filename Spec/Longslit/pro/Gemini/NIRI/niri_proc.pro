;+
; NAME:
;   niri_proc
;
; PURPOSE:
;   Overscan subtract (pattern noise remove) and flat field an NIRI 
;   spectroscopic image. 
;
; CALLING SEQUENCE:
;   niri_proc, filename, [ flux, invvar, hdr= $
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

pro niri_proc, filename, imag, invvar, hdr = hdr $
               , flatimg = flatimg, darkimg = darkimg $
               , adderr = adderr, verbose = verbose, bin = bin $
               , TELLURIC = TELLURIC, OSCAN = OSCAN

minval = 0.5
if (n_elements(adderr) EQ 0) then adderr = 0.01
hdr = xheadfits(filename)

if (size(hdr, /tname) NE 'STRING') then begin
    splog, 'Invalid FITS header for file ', filename
    imag = 0
    invvar = 0
    return
endif

;  If FLUX and INVVAR are not requested, then return with just the header
if (arg_present(imag) EQ 0 AND arg_present(invvar) EQ 0) then return

;Aladdin InSb (Hughes SBRC) Detector 
; See http://www.gemini.edu/sciops/instruments/niri/NIRIDetector.html

GAIN = 12.3D                     ; e-/ADU
RDNOISE_FAINT = 10.0D            ; e-/pi  (Low background mode) 
RDNOISE_TELL  = 70.0D            ; e-/pix (High background mode)

IF KEYWORD_SET(TELLURIC) THEN RDNOISE = RDNOISE_TELL $
ELSE RDNOISE = RDNOISE_FAINT

filelist = lookforgzip(filename)
if filelist[0] EQ '' then begin
    print, 'Could not find file named ', filename
    return
endif
hdr0 = xheadfits(filelist[0])
; lower wavelengths are down
imag = reverse(transpose(xmrdfits(filelist[0], 1)), 2)

if size(hdr0, /tname) EQ 'INT' then begin
    print, 'having trouble with ', filename
    return
endif

dims = size(imag, /dim)
nx = dims[0]
ny = dims[1]

imag = GAIN*imag
IF arg_present(invvar) THEN $
  invvar = 1.0/(abs(imag - sqrt(2.0)*rdnoise) +rdnoise^2)

IF KEYWORD_SET(OSCAN) THEN niri_oscan, imag, invvar = invvar

if (arg_present(invvar) AND adderr GT 0) then begin
    gmask = invvar GT 0         ; =1 for good points
    invvar = gmask / (1.0 / (invvar + (1-gmask)) + $
                      adderr^2 * (abs(imag))^2)
endif

ndim = size(imag, /n_dimen)
dims = size(imag, /dimens)

if (keyword_set(darkimg)) then begin
    if (size(darkimg, /n_dimen) NE ndim $
        OR total(size(darkimg, /dimens) NE dims) NE 0) then $
      message, 'Dimensions of image and darkimg image do not agree'
    
    imag = imag - darkimg
endif

if (keyword_set(flatimg)) then begin
    if (size(flatimg, /n_dimen) NE ndim $
        OR total(size(flatimg, /dimens) NE dims) NE 0) then $
      message, 'Dimensions of image and flat image do not agree'
    IF ARG_PRESENT(INVVAR) THEN invvar1 = invvar[20:976, *] 
    divideflat, imag[20:976, *], flatimg[20:976, *] $
                , invvar = invvar1, minval = minval $
                , quiet = (keyword_set(verbose) EQ 0)
;   For the moment, put the oscan back in
    
;   This will act like a bad pixel mask
    IF ARG_PRESENT(invvar) THEN BEGIN
        badmin = 0.7
        badmax = 1.3
        badpixmask = (flatimg[20:976, *] GT badmin AND flatimg LT badmax)
        invvar[20:976, *] = invvar[20:976, *]*badpixmask
    ENDIF
endif

return
end

;------------------------------------------------------------------------------
