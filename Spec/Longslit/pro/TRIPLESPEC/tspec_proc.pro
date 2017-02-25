;+
; NAME:
;   tspec_proc
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
;   tspec_clean  (from P. Muirhead, Cornell)
;
; REVISION HISTORY:
;   Apr-2011  Written by J. Hennawi (UCB), JXP
;-
;------------------------------------------------------------------------------
pro tspec_proc, filename, flux, invvar, hdr = hdr $
                , biasfile = biasfile, pixflatfile = pixflatfile $
                , illumflatfile = illumflatfile $
                , adderr = adderr, verbose = verbose, bin = bin $
                , rnoise = rnoise, gain = gain, maskbadpix = maskbadpix $
                , darkfile=darkfile
  
  if  N_params() LT 2  then begin 
     print, 'Syntax - ' + $
            'tspec_proc, fil, flux [invvar], HDR=, BIASFILE=, GAIN= [v1.0]' 
     return
  endif 



minval = 0.5
if (n_elements(adderr) EQ 0) then adderr = 0.01

hdr = xheadfits(filename)

if (size(hdr, /tname) NE 'STRING') then begin
stop
    splog, 'Invalid FITS header for file ', filename
    flux = 0
    invvar = 0
    return
endif

;  If FLUX and INVVAR are not requested, then return with just the header
if (arg_present(flux) EQ 0 AND arg_present(invvar) EQ 0) then return

;;  Clean image and convert to electrons
case strtrim(sxpar(hdr,'OBSERVAT'),2) of
   'APO': flux = tspec_apoclean(filename, hdr, ivar = invvar) 
   else: flux = tspec_clean(filename, hdr, ivar = invvar) 
endcase

;if (arg_present(invvar)) then $
;  tspec_oscan, filename, flux, invvar, hdr = hdr, verbose=verbose $
;              , rnoise = rnoise, gain = gain $
;else $
;   tspec_oscan, filename, flux, hdr = hdr, verbose = verbose $
;                , rnoise = rnoise, gain = gain

if (keyword_set(invvar) AND adderr GT 0) then begin
   gmask = invvar GT 0          ; =1 for good points
   invvar = gmask / (1.0 / (invvar + (1-gmask)) + $
                     adderr^2 * (abs(flux))^2)
endif
       
ndim = size(flux, /n_dimen)
dims = size(flux, /dimens)

;;;;;;;;;;;;;
;; Bias subtraction (electronics)

;; Removed bias-subtraction
;; This was the old unilluminated region, but it is now a dead amp
;; Took out bias subtraction 04/04/2012
;bias = djs_median(flux[7:37,*],1)
;savfilt = savgol(60,60,0,4 < (12 - 1))
;bias_fct = convol(bias, savfilt,/EDGE_TRUNCATE)
;bias_img = replicate(1., 1024L) # bias
;flux = flux - bias_img

;; This was code where I was playing around with creating a global
;; bias-correction using the unilluminated regions. 
;scatt_shift = 4L
;orderfile ='/Users/joe/REDUX_OTHER/P200_TSPEC/tspec-orders.fits'
;tset_slits = xmrdfits(orderfile, 1)
;xarr = findgen(dims[0]) # replicatE(1.0, dims[1])
;yarr = findgen(dims[1]) ## replicate(1.0, dims[0])
;slitmask = long_slits2mask(tset_slits)
;slitmask_left  = long_slits2mask(tset_slits, xshift = -abs(scatt_shift))
;slitmask_right = long_slits2mask(tset_slits, xshift = abs(scatt_shift))
;scattmask = (slitmask_left EQ 0 AND slitmask_right EQ 0 AND slitmask
;EQ 0)
;iscatt = where(scattmask)
;isort=sort(yarr[iscatt]) 
;ind = iscatt[isort]
;plot,yarr[ind],djs_median(flux[ind],width=1000),yr=[-10.0,30.0],/ysty,/xsty


;if (keyword_set(biasfile)) then begin
;    biasimg = xmrdfits(biasfile, 0, biashdr, $
;                      silent = (keyword_set(verbose) EQ 0))
;    if (size(biasimg, /n_dimen) NE ndim $
;        OR total(size(biasimg, /dimens) NE dims) NE 0) then $
;      message, 'Dimensions of image and bias image do not agree'
;    
;    flux = flux - biasimg
;endif

;;;;;;;;;;;;;
;; Dark subtraction (includes some K-band scattered light)
if (keyword_set(darkfile)) then begin
   dark = xmrdfits(darkfile,/silen)

   ;; Find scale values for dark frame
;   scale_sec = [1036, 1254, 858, 886]
   scale_sec = [858, 886, 791, 1012]
   sci_val = median(flux[scale_sec[0]:scale_sec[1], scale_sec[2]:scale_sec[3]])
   drk_val = median(dark[scale_sec[0]:scale_sec[1], scale_sec[2]:scale_sec[3]])
   scl_value = sci_val/drk_val

   ;; Gain kludge
;   g1 = median(flux[1000:1023,865:882] / flux[1024:1047,865:882]) 
   g1 = median(flux[865:882,1024:1047] / flux[865:882,1000:1023]) 
   flux[*,1024:*] = flux[*,1024:*] / g1 

   ;; Subtract off the dark (avoid K-band region)
   flux[0:894,*] = flux[0:894,*]-dark[0:894,*]*scl_value

   ;; A final kludge for the H-band?!

   ;; Check
;   xatv, flux, min=-50, max=500, /bloc
;   stop
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

;; Bad pixel mask
IF not (keyword_set(NO_MASKBADPIX)) then $
   tspec_badpixfix, flux, invvar

return
end

;------------------------------------------------------------------------------
