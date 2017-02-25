;+ 
; NAME:
; x_magtolum   
;    Version 1.1
;
; PURPOSE:
;  Convert AB mag to flux
;
; CALLING SEQUENCE:
;
; INPUTS:
;  inmag -- Magnitude
;  H0 -- Hubble constant
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
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Aug-2005 Written by JXP
;-
;------------------------------------------------------------------------------
function x_magtolum, inmag, h0

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'lum = x_magtolum(absmag, [h0]) (v1.0)'
      return, -1
  endif 

  ;; Convert mag according to h0
  if keyword_set(h0) then absmag = inmag - 5*alog10(h0) $
  else absmag = inmag

  ;; Flux
  flux = 10^(-1.d*(absmag + 48.6)/2.5)

  ;; Lum
  c = x_constants()
  lum = flux * 4 * !pi * (10*c.pc)^2
  stop

  return, lum
end


