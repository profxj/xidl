;+ 
; NAME:
; x_photsol1   
;   Version 1.1
;
; PURPOSE:
;    Performs the linear algebra on a set of obs with 1 free
;    parameters, i.e., Airmass fixed and NO color term
;
; CALLING SEQUENCE:
;   
; x_photsol1, mo0, mT, sig, airmass, setam, coeffs, sigma_coeff,
;              COLOR=
;
; INPUTS:
;   mo0 -  Observed magnitudes
;   mT  -  True magnitudes
;   sig - Combined Error 
;   airmass - Airmass values
;   SETAM - Value to set air mass coefficient to [Assumed positive!]
;
; RETURNS:
;   coeffs - Fit
;   sigma_coeffs - Error on the fit
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   color   - Color values
;
; OPTIONAL OUTPUTS:
;  CHISQ - chisq
;
; COMMENTS:
;
; EXAMPLES:
;   x_photsol1, blah
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_photsol1, mo0, mT, sig, airmass, setam, $
      coeffs, sigma_coeffs, CHISQ=chisq, COLOR=color

;
  if  N_params() LT 6  then begin 
      print, 'Syntax - ' +$
        'x_photsol1, mo0, mT, sig, airmass, setam, '
      print, '     coeffs, sigma_coeffs, CHISQ=, COLOR= (v 1.1)'
      return
  endif 

; Check setam
  if setam GT 0 then setam = -setam

; Optional keywords

  if keyword_set( COLOR ) then begin
      print, 'x_photsol1: Not ready for this yet!'
      return
  endif

  npoints = n_elements(mT)

; Weight is inverse variance 

   chk = where(sig EQ 0.0, cnt)
   if cnt GT 0 then begin
       print, 'Cant have zero error!'
       return
   endif

   weight        = 1.0/sig^2

;  Calculate weighted zeropoint

   zeropt = (mT-mo0-setam*airmass)*weight

   coeffs = total(zeropt)/total(weight)
   sigma_coeffs = sqrt(1./total(weight))
   
;  Chisq
   if arg_present( CHISQ ) then begin
       chisq = (total(((mT-mo0-setam*airmass) - coeffs)^2 * weight))/(npoints-1)
   endif

   return
end

