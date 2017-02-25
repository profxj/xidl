;+
; NAME:
;   long_proc
;
; PURPOSE:
;   Overscan subtract, bias subtract, and flat field a CCD image. The routine 
;   calls long_oscan and flat fields.
;
; CALLING SEQUENCE:
;   long_proc, filename, [ flux, invvar, hdr=, $
;    biasfile=, pixflatfile=, adderr=, /verbose ]
;
; INPUTS:
;   filename   -  name of a image file
;
; OPTIONAL INPUTS:
;   biasfile   - File with average bias
;   pixflatfile- File with pixel flat
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
; BUGS:
;   Need to add cosmic ray zapping.
;
; PROCEDURES CALLED:
;   divideflat
;   headfits()
;   mrdfits()
;   long_oscan
;
; REVISION HISTORY:
;   10-Mar-2005  Written by J. Hennawi (UCB), D. Schlegel (LBL)
;-
;------------------------------------------------------------------------------
pro mage_proc, filename, flux, invvar, hdr = hdr $
               , biasfile = biasfile, pixflatfile = pixflatfile $
               , illumflatfile=illumflatfile $
               , adderr = adderr, verbose = verbose, bin = bin $
               , rnoise = rnoise, gain = gain

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mage_proc, fil, flux [invvar], HDR=, BIASFILE=, GAIN= [v1.0]' 
      return
  endif 

minval = 0.5
if (n_elements(adderr) EQ 0) then adderr = 0.01
hdr = xheadfits(filename)

if (size(hdr, /tname) NE 'STRING') then begin
    splog, 'Invalid FITS header for file ', filename
    flux = 0
    invvar = 0
    return
endif

;  If FLUX and INVVAR are not requested, then return with just the header
if (arg_present(flux) EQ 0 AND arg_present(invvar) EQ 0) then return

if (arg_present(invvar)) then $
  mage_oscan, filename, flux, invvar, hdr = hdr, bin=bin, verbose=verbose $
              , rnoise = rnoise, gain = gain $
else $
   mage_oscan, filename, flux, hdr = hdr, verbose = verbose, bin = bin $
               , rnoise = rnoise, gain = gain

if (keyword_set(invvar) AND adderr GT 0) then begin
   gmask = invvar GT 0          ; =1 for good points
   invvar = gmask / (1.0 / (invvar + (1-gmask)) + $
                     adderr^2 * (abs(flux))^2)
endif
       
ndim = size(flux, /n_dimen)
dims = size(flux, /dimens)

if (keyword_set(biasfile)) then begin
    biasimg = xmrdfits(biasfile, 0, biashdr, $
                      silent = (keyword_set(verbose) EQ 0))
    if (size(biasimg, /n_dimen) NE ndim $
        OR total(size(biasimg, /dimens) NE dims) NE 0) then $
      message, 'Dimensions of image and bias image do not agree'
    
    flux = flux - biasimg
endif

if (keyword_set(pixflatfile)) then begin
    flatimg = xmrdfits(pixflatfile, 0, biashdr, $
                      silent = (keyword_set(verbose) EQ 0))
    if (size(flatimg, /n_dimen) NE ndim $
        OR total(size(flatimg, /dimens) NE dims) NE 0) then $
      message, 'Dimensions of image and flat image do not agree'
    
    divideflat, flux, flatimg, invvar = invvar, minval = minval, $
                quiet = (keyword_set(verbose) EQ 0)
 endif


if (keyword_set(illumflatfile)) then begin
   illumflat = xmrdfits(illumflatfile, 0, biashdr, $
                        silent = (keyword_set(verbose) EQ 0))
   if (size(illumflat, /n_dimen) NE ndim $
       OR total(size(illumflat, /dimens) NE dims) NE 0) then $
          message, 'Dimensions of image and illumflat image do not agree'
   ;; Don't apply illumination corrections larger than 30%
   gdpix = WHERE(illumflat GT 0.7D AND illumflat LT 1.3)
   flux[gdpix]  = flux[gdpix]/illumflat[gdpix]
   IF ARG_PRESENT(invvar) THEN $
      invvar[gdpix] = invvar[gdpix]*illumflat[gdpix]^2 
ENDIF

return
end

;------------------------------------------------------------------------------
