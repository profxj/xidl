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

pro nirspec_proc, filename, imag, invvar, hdr = hdr $
                  , pixflatfile =pixflatfile,darkfile=darkfile $
                  ,illumflatfile=illumflatfile $
                  , adderr = adderr, verbose = verbose, bin = bin
  
minval = 0.5
if (n_elements(adderr) EQ 0) then adderr = 0.01
hdr = xheadfits(filename,/silent)

if (size(hdr, /tname) NE 'STRING') then begin
    splog, 'Invalid FITS header for file ', filename
    imag = 0
    invvar = 0
    return
endif

;  If FLUX and INVVAR are not requested, then return with just the header
if (arg_present(imag) EQ 0 AND arg_present(invvar) EQ 0) then return

; Gain of NIRSEPC Aladdin-3 InSb, something seems wrong here, noise
; appears over-estimated. 
GAIN = 5.8                      ; e-/ADU
RDNOISE = 23.0                  ; e-
DARK_CURRENT = 0.8d             ; e-/s/pix
; For detector properties see:
; http://www2.keck.hawaii.edu/inst/nirspec/Specifications.html
exptime = float(sxpar(hdr, 'ITIME'))
dark_cnts = 0.8d*exptime
rdnoise_tot = sqrt(rdnoise^2 + dark_cnts)
;; Note we are not dealing with multiple co-adds etc. 

filelist = lookforgzip(filename)
if filelist[0] EQ '' then begin
   print, 'Could not find file named ', filename
   return
endif
hdr0 = xheadfits(filelist[0],/silent)
imag = xmrdfits(filelist[0], 0, /silent)

if size(hdr0, /tname) EQ 'INT' then begin
    print, 'having trouble with ', filename
    return
endif

dims = size(imag, /dim)
nx = dims[0]
ny = dims[1]

imag = GAIN*imag
IF arg_present(invvar) THEN $
   invvar = 1.0/(abs(imag - sqrt(2.0)*rdnoise_tot) +rdnoise_tot^2)

if (arg_present(invvar) AND adderr GT 0) then begin
    gmask = invvar GT 0         ; =1 for good points
    invvar = gmask / (1.0 / (invvar + (1-gmask)) + $
                      adderr^2 * (abs(imag))^2)
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
