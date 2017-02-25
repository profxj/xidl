;+
; NAME:
;   not_proc
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
pro not_proc, filename, flux, invvar, ihdr = ihdr, hdr_arr = hdr_arr $
              , superbias = superbias, delta_bias = delta_bias $
              , superflat = superflat, instbiasfile = instbiasfile $
              , superdark = superdark $
              , adderr = adderr, verbose = verbose, bin = bin $
              , rnoise = rnoise, gain = gain, mask = mask

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'long_proc, fil, flux [invvar], HDR=, BIASFILE=, GAIN= [v1.0]' 
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

not_biassub, filename, flux, invvar, superbias = superbias $
             , delta_bias = delta_bias, instbiasfile = instbiasfile $
             , superdark = superdark $
             , hdr_arr = hdr_arr, mask = mask, ihdr = ihdr

if (keyword_set(invvar) AND adderr GT 0) then begin
   gmask = invvar GT 0          ; =1 for good points
   invvar = gmask / (1.0 / (invvar + (1-gmask)) + $
                     adderr^2 * (abs(flux))^2)
endif
       
ndim = size(flux, /n_dimen)
dims = size(flux, /dimens)

IF KEYWORD_SET(superflat) THEN BEGIN
   if (size(superflat, /n_dimen) NE ndim $
       OR total(size(superflat, /dimens) NE dims) NE 0) then $
          message, 'Dimensions of image and flat image do not agree'
   divideflat, flux, superflat, invvar = invvar, minval = minval, $
               quiet = (keyword_set(verbose) EQ 0)
endif

return
end

;------------------------------------------------------------------------------
