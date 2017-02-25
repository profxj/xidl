;+ 
; NAME:
; x_photsol2   
;   Version 1.1
;
; PURPOSE:
;    Performs the linear algebra on a set of obs with 2 free
;    parameters, i.e., either Airmass fixed or NO color term
;
; CALLING SEQUENCE:
;   
; x_photsol2, mo0, mT, sig, airmass, coeffs, sigma_coeff,
;              SETAM=, COLOR=, CHISQ=
;
; INPUTS:
;   mo0 -  Observed magnitudes
;   mT  -  True magnitudes
;   sig - Combined Error 
;   airmass - Airmass values
;
; RETURNS:
;   coeffs - Fit
;   sigma_coeffs - Error on the fit
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   SETAM - Value to set air mass coefficient to [Assumed positive!]
;   color   - Color values
;
; OPTIONAL OUTPUTS:
;  CHISQ - chisq
;
; COMMENTS:
;
; EXAMPLES:
;   x_photsol2, blah
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_photsol2, mo0, mT, sig, airmass, $
      coeffs, sigma_coeffs, CHISQ=chisq, SETAM=setam, COLOR=color

;
  if  N_params() LT 6  then begin 
      print, 'Syntax - ' +$
        'x_photsol2, mo0, mT, sig, airmass, '
      print, '     coeffs, sigma_coeffs, CHISQ=, SETAM=, COLOR= (v 1.0)'
      return
  endif 

; Optional keywords

  if keyword_set( SETAM ) then begin
      if not keyword_set( COLOR ) then begin
          print, 'x_photsol2: Must provide color when AM is set!'
          return
      endif
      if setam LT 0 then begin
          print, 'x_photsol2: SETAM assumed positive!', setam
          return
      endif
  endif

;  Assumes 2 terms: 

   nterms = 2
   npoints = n_elements(mo0)
 
; weight is inverse variance 

   chk = where(sig EQ 0.0, cnt)
   if cnt GT 0 then begin
       print, 'Cant have zero error!'
       return
   endif

   weight        = 1.0/sig^2

   action_matrix = fltarr(npoints,nterms) + 1.0

   
; Set action to Airmass or Color term as appropriate

   if keyword_set( SETAM ) then begin
       ; Offest the observed magnitudes by the airmass term!
       mo0 = mo0 - setam*airmass
       action_matrix[*,1] = -1.*color
   endif else action_matrix[*,1] = -1.*airmass



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

   return
end

