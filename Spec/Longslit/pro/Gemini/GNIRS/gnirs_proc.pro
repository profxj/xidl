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

pro gnirs_proc, filename, imag, invvar, hdr = hdr $
                , flatimg = flatimg, darkimg = darkimg $
                , adderr = adderr, verbose = verbose, bin = bin $
                , TELLURIC = TELLURIC, PATTERN = PATTERN

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

; Gain of GNIRS Aladdin II detector
GAIN = 13.0  ; e-/ADU
RDNOISE_FAINT = 7.0             ; e-
RDNOISE_TELL = 38.0             ; e-
; For detector properties see:
; http://gemini.ast.cam.ac.uk//sciops/instruments/nirs/nirsDetector.html

IF KEYWORD_SET(TELLURIC) THEN RDNOISE = RDNOISE_TELL $
ELSE RDNOISE = RDNOISE_FAINT

filelist = lookforgzip(filename)
if filelist[0] EQ '' then begin
    print, 'Could not find file named ', filename
    return
endif
hdr0 = xheadfits(filelist[0])
; lower wavelengths are down
imag = xmrdfits(filelist[0], 1, /silent)

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

if (arg_present(invvar) AND adderr GT 0) then begin
    gmask = invvar GT 0         ; =1 for good points
    invvar = gmask / (1.0 / (invvar + (1-gmask)) + $
                      adderr^2 * (abs(imag))^2)
;;    invvar[*, 1010:*] = 0.0D    ; top of chip is not illuminated
endif

ndim = size(imag, /n_dimen)
dims = size(imag, /dimens)

if (keyword_set(flatimg)) then begin
    if (size(flatimg, /n_dimen) NE ndim $
        OR total(size(flatimg, /dimens) NE dims) NE 0) then $
      message, 'Dimensions of image and flat image do not agree'
;   Don't mask the flat
    divideflat, imag, flatimg, invvar = invvar, minval = 0.0 $
                , quiet = (keyword_set(verbose) EQ 0)
ENDIF

IF KEYWORD_SET(PATTERN) THEN gnirs_remove_patternnoise, imag

return
end

;------------------------------------------------------------------------------
