;+ 
; NAME:
; mike_gain
;     Version 2.0
;
; PURPOSE:
;    Calculate the inverse gain, e-/DN with two very similiar input images
;      For MIKE, milky flats or well illuminated trace flats are best.
;
; CALLING SEQUENCE:
;  estimated_gain = mike_gain(image1, image2)
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
;      gain = mike_gain(image1, image2)
;
; PROCEDURES/FUNCTIONS CALLED:
;   x_calcgain
;
; REVISION HISTORY:
;   16-Jul-2003 Written by SMB
;   20-Jun-2005 Consumed into XIDL by JXP
;-
;------------------------------------------------------------------------------

function mike_gain, image1, image2

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_gain(image1, image2) (v2.0)'
      return, -1
  endif 

  ;; x_calcgain
  gain = x_calcgain(image1, image2)
  return, gain

end

