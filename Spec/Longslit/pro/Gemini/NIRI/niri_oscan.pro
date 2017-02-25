; NAME:
;   niri_oscan
;
; PURPOSE:
;   Overscan subtracts NIRI spectroscopic data removing pattern noise
;
; CALLING SEQUENCE:
;   niri_oscan, filename, [ rawsub, rawivar, hdr= , $
;    gain= , rnoise= , image=, /verbose ]
;
; INPUTS:
;   filename   -  name of a image file
;
; OPTIONAL INPUTS:
;   verbose    - If set, then verbose output
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   rawsub     - Bias subtracted image
;   rawivar    - Inverse variance of image
;   hdr        - Fits header of raw image
;   gain       - vector of gains for each amp
;   rnoise     - vector or readnoise for each amps
;   image      - original data image
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   mrdfits
;   djs_avsigclip
; REVISION HISTORY:
;   10-Mar-2005  Written by J. Hennawi (UCB)
;-
;------------------------------------------------------------------------------


pro niri_oscan, imag, invvar = invvar, inmask = inmask1

; inmask has avsigclip convention
oscan_iml = imag[0:10, *]
IF KEYWORD_SET(INMASK1) THEN inmask = inmask1[0:10, *]
oscan_l = djs_avsigclip(oscan_iml, 1, inmask = inmask, sigrej = 2.0)
oscan_imr = imag[981:1020, *]
IF KEYWORD_SET(INMASK1) THEN inmask  = inmask1[981:1020, *]
oscan_r = djs_avsigclip(oscan_imr, 1, inmask = inmask, sigrej = 2.0)

imag[0:511, *] = imag[0:511, *] - replicate(1.0, 512) # oscan_l
imag[512:*, *] = imag[512:*, *] - replicate(1.0, 512) # oscan_r
; Zero the overscan region in invvar imag
IF KEYWORD_SET(invvar) THEN BEGIN
    invvar[0:19, *] = 0.0
    invvar[976:*, *] = 0.0
ENDIF

RETURN
END


