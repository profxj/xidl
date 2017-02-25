;+ 
; NAME:
; x_fluxstar
;   Version 1.0
;
; PURPOSE:
;    Calcualtes the flux at surface of a star assuming Plank's
;    law for an energy or energy interval
;
; CALLING SEQUENCE:
;   
;  flux = x_fluxstar(energy, T)
;
; INPUTS:
;  energy -  Value or range (2 values) in eV
;  T      -  Temperature of the star
;
; RETURNS:
;   
; OUTPUTS:
;   flux (ergs/s/cm^2)  [Hz^-1 for single energy]
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
;   12-Dec-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_fluxstar, energy, T

  if  N_params() LT 2  then begin 
      print, 'Syntax - ' +$
        'flux = x_fluxstar( energy {eV}, T )  (v1.1)'
      return, -1
  endif 

  nstp = 10000L
  c = x_constants()

  ;;
  if n_elements(energy) EQ 1 then begin
      ;; Planck
      lambda = c.h * c.c / (energy * c.eV)  ; cm
      lambda = lambda * 1e8  ; Ang
      flux = planck(lambda, T)  ; ergs/s/cm^2/Ang
  endif else begin
      ;; dLog E
      deng = alog(energy[1]/energy[0]) / (nstp-1.)
      engs = exp( alog(energy[0]) + dindgen(nstp)*deng)
      ;; Planck
      lambda = c.h * c.c / (engs * c.eV)  ; cm
      lambda = lambda * 1e8 ; Ang
      fluxs = planck(lambda, T)  ; F_lambda
      nus = engs * c.eV / c.h  ; nu
      lnus = alog(nus)
      dlnus = lnus - shift(lnus,1) ;  d ln nu
      dlnus[0] = dlnus[1]
      ;; Total
      fluxnu = fluxs * lambda^2 / c.c / 1e8  ;  Fnu
      flux = total( nus * fluxnu * dlnus )  ; ergs/s/cm^2
  endelse

return, flux
end

