;+ 
; NAME:
; x_photsol3   
;   Version 1.1
;
; PURPOSE:
;    Performs the linear algebra on a set of obs.  Fits for three
;  free parameters: zero point, airmass term, color term
;
; CALLING SEQUENCE:
;   
; x_photsol3, mo0, mT, sig, airmass, color, coeffs, sigma, 
;      CHISQ=, NCORR=
;
; INPUTS:
;   mo0 -  Observed magnitudes
;   mT  -  True magnitudes
;   sig - Combined Error 
;   airmass - Airmass values
;   color   - Color values
;
; RETURNS:
;   coeffs - Fit
;   sigma_coeffs - Error on the fit
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  NCORR - Normalized correlation matrix
;  CHISQ - chisq
;
; COMMENTS:
;
; EXAMPLES:
;   x_photsol3, blah
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_photsol3, mo0, mT, sig, airmass, color, $
      coeffs, sigma_coeffs, CHISQ=chisq, NCORR=ncorr

;
  if  N_params() LT 7  then begin 
      print, 'Syntax - ' +$
        'x_photsol3, mo0, mT, sig, airmass, color, '
      print, '     coeffs, sigma_coeffs, CHISQ=, NCORR='
      return
  endif 

;  Assumes 3 terms: 

   nterms = 3
   npoints = n_elements(mo0)
 
; weight is inverse variance 

   chk = where(sig EQ 0.0, cnt)
   if cnt GT 0 then begin
       print, 'Cant have zero error!'
       return
   endif

   weight        = 1.0/sig^2

   action_matrix = fltarr(npoints,nterms) + 1.0

   ; 0th row is zero_point, represented by all 1's
   ; 1st row is airmass term, 2nd is color...

   action_matrix[*,1] = -1.*airmass
   action_matrix[*,2] = -1.*color

   alpha = transpose(action_matrix) # $
         (action_matrix * (weight # replicate(1,nterms)))
   beta  = transpose(((mT-mo0) * weight) # action_matrix)

   correlation_matrix = invert(alpha, status, /double)
   if status NE 0 then begin
       print, 'x_photsol: Probably a singular matrix', status
       return
   endif

   coeffs = correlation_matrix # beta
   sigma_coeffs = sqrt(correlation_matrix[lindgen(nterms),lindgen(nterms)])
   
;  Chiqs
   if arg_present( CHISQ ) then begin
       modelfit = action_matrix # coeffs
       chisq = (total(((mT-mo0) - modelfit)^2 * weight))/(npoints-nterms)
   endif

;  Normalized correlation matrix
   if keyword_set( NCORR ) then begin
       norm = sigma_coeffs ## sigma_coeffs
       ncorr = correlation_matrix / norm
   endif
   
   return
end

