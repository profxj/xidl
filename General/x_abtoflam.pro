;+ 
; NAME:
; x_abtoflam
;
; PURPOSE:
;    Converts AB mag to flambda (given a wavelength)
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;  AB   -- AB magnitude
;  WAVE -- Wavelength (Angstroms)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   24-Jan-2008 Written by JXP
;-
;------------------------------------------------------------------------------
function x_abtoflam, AB, wave, FNU=fnu

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'x_abtoflam, AB, wave, FNU= (v1.0)'
    return, -1
  endif 

  ;; Convert to FNU
  fnu = 10.d^( -1. * (AB+48.6) * 0.4 )  ;; erg s^-1 cm^-2 Hz^-1

  ;; Convert to flambda
  c = x_constants()
  flam = fnu * c.c / wave / (wave * 1e-8)  ;; erg s^-1 cm^-2 A^-1

  return, flam

end

