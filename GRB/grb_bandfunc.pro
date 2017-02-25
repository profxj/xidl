;+ 
; NAME:
;  grb_fxlum
;   Version 1.1
;
; PURPOSE:
;    Calculate the flux for a Band function.
;
; CALLING SEQUENCE:
;   
;   grb_fxlum, fini, tini, freq, RAD=, FLUX=, LUM=, Z=, $
;              tEarth=,  TEXP=, TSTART=, NSTP=
;
; INPUTS:
;       E -- Energy (keV)
;   alpha -- Power-law exponent
;   beta  -- Power-law exponent
;
; RETURNS:
;   flux= -- Array of flux values computed as function of time
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;    A=     -- Radiative expansion [Default: adiabatic]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   21-Apr-2008 Written by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function grb_bandfunc, E, Epeak, alpha, beta, A=A 

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'f(E) = grb_bandfunc(E,Epeak, [alpha, beta], A=A) [v1.1]'
    return, -1
  endif 

  ;; Optional keywords
  if not keyword_set( alpha ) then alpha = -1.
  if not keyword_set( beta ) then beta = -2.
  if not keyword_set( A ) then a = 1.  ;; photons / s / cm^2 / keV

  fE = fltarr(n_elements(E))  ;; Photon flux

  c = x_constants()

  ;; Functional form
  Ec = (alpha-beta)*Epeak/(2+alpha)  ;; keV

  low = where(E LT Ec, nlow, complement=hi, ncomplement=nhi)

  ;; Low
  if nlow NE 0 then fE[low] = A*(E[low]/100.)^alpha * $
    exp(-E[low]*(2+alpha)/Epeak)
  if nhi NE 0 then fE[hi] = A* (Ec/100)^(alpha-beta) * $
    exp(beta-alpha) * (E[hi]/100)^beta


  return, fE
end
