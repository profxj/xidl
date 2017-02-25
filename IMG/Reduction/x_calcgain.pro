;+ 
; NAME:
; x_calcgain
;     Version 1.1
;
; PURPOSE:
;    Calculate the inverse gain, e-/DN with two very similiar input images
;      For MIKE, milky flats or well illuminated flats are best.
;
; CALLING SEQUENCE:
;   
;  estimated_gain = x_gain(image1, image2)
;
; INPUTS:
;   image1  -  Any array of pixels with reasonable gaussian errors
;   image2  -  A very similar image to image1, but not identical!
;
; RETURNS:
;  gain   -  A floating point number: the inverse gain
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  MINPIX -- Minimum number of pixels with DN>MINVAL in both images to
;            perform statistics  (Default: 1000)
;  MINVAL -- Minimum flux to perform stats (default: 100)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The gain is calculated from the variance of the natural logarithm 
;     of the ratio of two images which have units of ADU (Poisson counts / gain)
;
;     <gain> = 1.0 / { variance [ sqrt(reduced_flux) * ln(i1/i2) ] } 
;
; EXAMPLES:
;      image1 = readfits('ov_mr0405.fits') 
;      image2 = readfits('ov_mr0406.fits')
;      gain = x_calcgain(image1, image2)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   16-Jul-2003 Written by SMB
;   04-Apr-2005 Ported to XIDL by JXP
;-
;------------------------------------------------------------------------------

function x_calcgain, image1, image2, MINPIX=minpix, MINVAL=minval

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'x_calcgain(image1, image2, MINPIX=, MINVAL=minval) (v2.0)'
      return, -1
  endif 
  
  if not keyword_set(MINPIX) then minpix = 1000L
  if not keyword_set(MINVAL) then minval = 100.

  ;;
  i1 = transpose(image1)
  i2 = transpose(image2)
  
  ;; Grab the good pixels
  pos = where(i1 GT MINVAL AND i2 GT MINVAL, npos)
  if npos LT MINPIX then begin
      print, 'x_calcgain: Not enoung positive pixel pairs'
      return, 0.0
  endif
  
  ;; Grab the ratio
  x =  i1[pos]/i2[pos]
  r = median(x)
  print, 'x_calcgain: Ratio of images is ', r
  reduced_flux = 1.0 / (1./i1[pos]  + 1./i2[pos])
  
  ;; Guess
  guess = sqrt(reduced_flux)* alog(x/r)
  sig_guess = stddev(guess)
  keep = where(abs(guess) LT 4.0*sig_guess AND $
               reduced_flux GT median(reduced_flux) , nkeep)
  
  if nkeep LT 100 then begin
      print, 'x_calcgain: Too many outliers'
      return, 0.0
  endif
  
  
  ;; Smooth
  smooth_guess = smooth(guess[keep],11) 
  
  ;; What is the 11/10 for??  JXP
  gain = 1./(mean((guess[keep]-smooth_guess)^2)*11./10.)
  print, 'x_calcgain: The gain is ', gain


  return, gain
end

