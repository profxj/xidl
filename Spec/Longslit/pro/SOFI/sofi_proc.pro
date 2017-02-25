;+
; NAME:
;   sofi_proc
;
; PURPOSE:
;   Process a SOFI near-IR image
;
; CALLING SEQUENCE:
;   sofi_proc, filename, [ flux, invvar, hdr= $
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

pro sofi_proc, filename, imag, invvar, hdr = hdr $
               , pixflatfile = pixflatfile, darkfile = darkfile $
               , illumflatfile = illumflatfile $
               , adderr = adderr, verbose = verbose, bin = bin

 
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

;; SOFI Rockwell Hg:Cd:Te 1024x1024 Hawaii a
GAIN = 5.4 ;; e-/ADU
NDIT = float(esopar(hdr, 'HIERARCH ESO DET NDIT '))
IF NDIT LE 0.0 THEN NDIT = 1.0
ind_samp = WHERE(stregex(hdr, 'HIERARCH ESO DET NCORRS NAME', /bool))
sampling = strcompress(strmid(hdr[ind_samp], 30, 14))
IF strmatch(sampling, '*Fowler*') THEN BEGIN
   ind_nsamp = WHERE(stregex(hdr, 'HIERARCH ESO DET NDSAMPLES', /bool))
   nsamp = (double(strmid(hdr[ind_nsamp], 30, 14)))[0]
   CASE nsamp OF 
      4:  RDNOISE = 7.63
      8:  RDNOISE = 6.40
      10: RDNOISE = 5.81
      12: RDNOISE = 5.24
      16: RDNOISE = 4.88
      20: RDNOISE = 4.52
      30: RDNOISE = 4.23
      ELSE: BEGIN 
         ;; !! Bug for NSAMP < 4 or nsamp > 30 it will extrapolate
         samp_vec = [4.00, 8.00, 10.0, 12.0, 16.0, 20.0, 30.0]
         rn_vec =   [7.63, 6.40, 5.81, 5.24, 4.88, 4.52, 4.23]
         RDNOISE = interpol(rn_vec, samp_vec, float(nsamp))
      END 
   ENDCASE
ENDIF ELSE RDNOISE = 10.0d ;; measure this
filelist = lookforgzip(filename)
if filelist[0] EQ '' then begin
    print, 'Could not find file named ', filename
    return
endif
hdr0 = xheadfits(filelist[0])
; lower wavelengths are up by default, so reverse to flip
imag = reverse(xmrdfits(filelist[0], 0, /silent), 2)


if size(hdr0, /tname) EQ 'INT' then begin
    print, 'Having trouble with ', filename
    return
endif

dims = size(imag, /dim)
nx = dims[0]
ny = dims[1]

imag = GAIN*NDIT*imag
IF arg_present(invvar) THEN $
   invvar = 1.0/(abs(imag - sqrt(2.0)*rdnoise) +rdnoise^2)
;; Add bad pixel mask here
if (arg_present(invvar) AND adderr GT 0) then begin
   gmask = invvar GT 0          ; =1 for good points
   invvar = gmask / (1.0 / (invvar + (1-gmask)) + $
                     adderr^2 * (abs(imag))^2)
endif

ndim = size(imag, /n_dimen)
dims = size(imag, /dimens)

if (keyword_set(darkfile)) then begin
   darkhdr = xheadfits(darkfile)
   stacked_dark = long(sxpar(darkhdr, 'STACKED'))
   IF stacked_dark THEN $
      darkimg = xmrdfits(darkfile, 0, darkhdr, $
                         silent = (keyword_set(verbose) EQ 0)) $
   ELSE sofi_proc, darkfile, darkimg, hdr = darkhdr
   if (size(darkimg, /n_dimen) NE ndim $
       OR total(size(darkimg, /dimens) NE dims) NE 0) then $
          message, 'Dimensions of image and darkimg image do not agree'
   imag = imag - darkimg
endif

if (keyword_set(pixflatfile)) then begin
   flatimg = xmrdfits(pixflatfile, 0, flathdr, $
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
   illumflat = xmrdfits(illumflatfile, 0, illumhdr, $
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
