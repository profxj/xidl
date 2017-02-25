;+ 
; NAME:
; x_coolrate
;   Version 1.0
;
; PURPOSE:
;    Calculates the cooling rate for a T<10000K gas assuming Lya and
;    fine-structure lines.  Uses equations from Hollenbach & McKee
;    (1989)
;
; CALLING SEQUENCE:
;   
; cool_rate = x_coolrate( nH, x, T, [Z])
;
; INPUTS:
;   nH     - Volume density of H
;   x      - Ionization fraction  (n_e = x * n_H)
;   T      - Gas temperature
;   [Z]      - Metallicity
;
; RETURNS: 
;  cool_rate -- Cooling rate
;   
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Dec-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_coolrate, nH, x, T, Z

  if  N_params() LT 3  then begin 
      print, 'Syntax - ' +$
        'cool_rate = x_coolrate( nH, x, T, Z ) (v1.1)'
      return, -1
  endif 
;
  ;; Optional keywords
  if not keyword_set(Z) then Z = 1.
  c = x_constants()

  nT = n_elements(T)
  cool_rate = dblarr(nT)

  ;; Line cooling
  parse_atmcool, atom_cool, Z=Z

  for qq=0L,nT-1 do begin
      T2 = T[qq] / 100.d
      ;; Lya
      LambdaLya = 7.3d-19 * x * exp(-1184./T2)
      
      
      gd = where(atom_cool.Ti GT 0.)
      
      gd_atm = atom_cool[gd]
      
      ;; Values
      degen = (2* gd_atm.jh + 1) / (2*gd_atm.jl + 1.)
      dencrite = gd_atm.ae * T2^gd_atm.xe
      dencrith = gd_atm.ah * T2^gd_atm.xh
      
      anume = gd_atm.gae * T2^gd_atm.gxe * degen * exp(-1*gd_atm.Ti/T[qq])
      anumh = gd_atm.gah * T2^gd_atm.gxh * degen * exp(-1*gd_atm.Ti/T[qq])
      adenom = 1. + nH * (x/dencrite + (1.-x)/dencrith)
      
      ALe = x * c.k * gd_atm.Ti * anume / adenom
      ALh = (1.-x) * c.k * gd_atm.Ti * anumh / adenom

      Alambda = total( gd_atm.abnd * (ALe + ALh) )
      
      ;; Save
      cool_rate[qq] = nH * (Alambda + LambdaLya)
  endfor
      
  return, cool_rate
end
