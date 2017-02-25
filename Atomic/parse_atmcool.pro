;+ 
; NAME:
; parse_atmcool
;   Version 1.0
;
; PURPOSE:
;    Parse the atomic data for cooling
;
; CALLING SEQUENCE:
;   parse_atmcool, struct, Z=Z
;
; INPUTS:
;
; RETURNS: 
;   
; OUTPUTS:
;  atom_cool -- Structure holding the cooling rates
;
; OPTIONAL OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  Z=  -- Metallicity of the gas [Default = 1]
;
; COMMENTS:
;
; EXAMPLES:
;   
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Dec-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro parse_atmcool, atom_cool, Z=Z

  if  N_params() LT 1  then begin 
      print, 'Syntax - ' +$
        'parse_atmcool, atom_cool, Z= [v1.0]'
      return
  endif 
;
  ;; Optional keywords
  if not keyword_set(Z) then Z = 1.

  tmp = { $
          lbl: '', $
          abnd: 0., $  ; Abundance
          Ti: 0., $             ; Energy (K)
          ae: 0., $             ; Coefficient for electron critical density
          xe: 0., $             ; Exponent for electron critical density
          ah: 0., $             ; Coefficient for electron critical density
          xh: 0., $             ; Exponent for electron critical density
          gae: 0., $            ; Coefficient for electron recombination
          gxe: 0., $            ; Exponent for electron recombination
          gah: 0., $            ; Coefficient for hydrogen recombination
          gxh: 0., $            ; Exponent for hydrogen recombination
          jl: 0., $             ; j value of lower state
          jh: 0. $              ; j value of higher state
        }
  atom_cool = replicate(tmp, 10)    ; 0=CII,

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Carbon
  readcol, getenv('XIDL_DIR')+'/Atomic/Data/CIIcool.dat', $
    lbl, CII_Ti, CII_ae, CII_xe, CII_ah, CII_xh, CII_gae, CII_gxe, $
    CII_gah, CII_gxh, CII_jl, CII_jh, FORMAT='A,F,F,F,F,F,F,F,F,F,F'

  getabnd, 'C', ZZ, abnd

  atom_cool[0].lbl = lbl
  atom_cool[0].abnd = Z* (10.^(abnd-12.))
  atom_cool[0].Ti = CII_Ti
  atom_cool[0].ae = CII_ae
  atom_cool[0].xe = CII_xe
  atom_cool[0].ah = CII_ah
  atom_cool[0].xh = CII_xh
  atom_cool[0].gae = CII_gae
  atom_cool[0].gxe = CII_gxe
  atom_cool[0].gah = CII_gah
  atom_cool[0].gxh = CII_gxh
  atom_cool[0].jl = CII_jl
  atom_cool[0].jh = CII_jh

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Iron
  readcol, getenv('XIDL_DIR')+'/Atomic/Data/FeIIcool.dat', $
    lbl, FeII_Ti, FeII_ae, FeII_xe, FeII_ah, FeII_xh, FeII_gae, FeII_gxe, $
    FeII_gah, FeII_gxh, FeII_jl, FeII_jh, FORMAT='A,F,F,F,F,F,F,F,F,F,F'

  getabnd, 'Fe', ZZ, abnd

  atom_cool[1:2].lbl = 'FeII'
  atom_cool[1:2].abnd = Z* (10.^(abnd-12.))
  atom_cool[1:2].Ti = FeII_Ti
  atom_cool[1:2].ae = FeII_ae
  atom_cool[1:2].xe = FeII_xe
  atom_cool[1:2].ah = FeII_ah
  atom_cool[1:2].xh = FeII_xh
  atom_cool[1:2].gae = FeII_gae
  atom_cool[1:2].gxe = FeII_gxe
  atom_cool[1:2].gah = FeII_gah
  atom_cool[1:2].gxh = FeII_gxh
  atom_cool[1:2].jl = FeII_jl
  atom_cool[1:2].jh = FeII_jh

  return
end
