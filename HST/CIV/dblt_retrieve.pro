;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; dblt_find.pro                
; Author: Kathy Cooksey                      Date: 4 Nov 2005
; Project: 
; Description: 
;   Rest wavelengths and oscillator strengths from lls.lst
;   and gamma values from atom.dat.
; Input: 
; Optional:
; Output: 
; Example:
; History:
;   4 Nov 2005  Created KLC
;  29 May 2011  Update
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function dblt_retrieve,dblt_name

  neviii = { doubletstrct }
  neviii.ion = 'NeVIII'
  neviii.wvI = 770.409
  neviii.fI = 0.1030
  neviii.gammaI = 0.d
  neviii.wvII = 780.324
  neviii.fII = 0.0505
  neviii.gammaII = 0.d

  ovi = { doubletstrct }
  ovi.ion = 'OVI'
  ovi.wvI = 1031.9261
  ovi.fI = 0.13290
  ovi.gammaI = 4.163e8 
  ovi.wvII = 1037.6167
  ovi.fII = 0.06609 
  ovi.gammaII = 4.095e8

  nv = { doubletstrct }
  nv.ion = 'NV'
  nv.wvI = 1238.821
  nv.fI = 0.1570
  nv.gammaI = 3.411e8 
  nv.wvII = 1242.804
  nv.fII = 0.078230 
  nv.gammaII = 3.378e8

  siiv = { doubletstrct }
  siiv.ion = 'SiIV'
  siiv.wvI = 1393.755 
  siiv.fI = 0.5280
  siiv.gammaI = 8.825e8
  siiv.wvII = 1402.770
  siiv.fII = 0.262 
  siiv.gammaII = 8.656e8

  fuiv = { doubletstrct }
  fuiv.ion = 'FuIV'             ; faux doublet for testing
  fuiv.wvI = 1393.755 - 0.6
  fuiv.fI = 0.5280
  fuiv.gammaI = 8.825e8
  fuiv.wvII = 1402.770 + 1.7
  fuiv.fII = 0.262 
  fuiv.gammaII = 8.656e8

  civ = { doubletstrct }
  civ.ion = 'CIV'
  civ.wvI = 1548.195
  civ.fI = 0.19080
  civ.gammaI = 2.654E8
  civ.wvII = 1550.770
  civ.fII = 0.095220 
  civ.gammaII = 2.641E8 

  fkiv = { doubletstrct }
  fkiv.ion = 'FkIV'             ; fake doublet for testing
  fkiv.wvI = 1548.195 - 0.6
  fkiv.fI = 0.19080
  fkiv.gammaI = 2.654E8
  fkiv.wvII = 1550.770 + 1.7
  fkiv.fII = 0.095220 
  fkiv.gammaII = 2.641E8 

  lya = { doubletstrct }
  lya.ion = 'Lya' 
  lya.wvI = 1215.6701
  lya.fI = 0.41640 
  lya.gammaI = 6.265E8
  lya.wvII = 1025.7223
  lya.fII = 0.079120
  lya.gammaII = 1.897E8  

  ciii = { doubletstrct }
  ciii.ion = 'CIII'
  ciii.wvI = 977.020
  ciii.fI = 0.7620
  ciii.gammaI = 1.775E9
  ciii.wvII = lya.wvI
  ciii.fII = lya.fI
  ciii.gammaII = lya.gammaI

  caii = { doubletstrct }
  caii.ion = 'CaII'
  caii.wvI = 3934.777
  caii.fI = 0.650
  caii.gammaI = 1.456E8
  caii.wvII = 3969.591
  caii.fII = 0.322
  caii.gammaII = 1.414E8

  mgii = { doubletstrct }
  mgii.ion = 'MgII'
  mgii.wvI = 2796.352
  mgii.fI = 0.6123
  mgii.gammaI = 2.612E8
  mgii.wvII = 2803.531
  mgii.fII = 0.3054
  mgii.gammaII = 2.592E8

  fkii = { doubletstrct }
  fkii.ion = 'FkII'             ; Slight offset from MgII sep.
  fkii.wvI = 2796.352 - 0.6
  fkii.fI = 0.6123
  fkii.gammaI = 2.612E8
  fkii.wvII = 2803.531 - 1.7
  fkii.fII = 0.3054
  fkii.gammaII = 2.592E8

  feii = { doubletstrct }
  feii.ion = 'FeII'
  feii.wvI = 2600.1729
  feii.fI = 0.2390
  feii.gammaI = 2.700E8
  feii.wvII = 2586.650
  feii.fII = 0.0691
  feii.gammaII = 2.720E8

  if size(dblt_name,/type) eq 7 then begin
     name = strlowcase(strtrim(dblt_name,2))
     case name of 
        'neviii': return,neviii
        'ovi': return,ovi
        'nv': return,nv
        'siiv': return,siiv
        'fuiv': return,fuiv
        'civ': return,civ
        'fkiv': return,fkiv     ; fake
        'lya': return,lya
        'hi': begin
           lya.ion = 'HI' 
           return, lya
        end
        'ciii': return, ciii
        'caii': return, caii
        'mgii': return, mgii
        'fkii': return, fkii
        'feii': return, feii
        else: return,{ doubletstrct }
     endcase
  endif else begin
     ;; Assume it's a number and have to find it
     info = [ $
            ['neviii',string(neviii.wvi)],$
            ['ovi',string(ovi.wvi)],$
            ['nv',string(nv.wvi)],$
            ['siiv',string(siiv.wvi)],$
            ['fuiv',string(fuiv.wvi)],$
            ['civ',string(civ.wvi)],$
            ['fkiv',string(fkiv.wvi)],$
            ['lya',string(lya.wvi)],$
            ['hi',string(lya.wvi)],$
            ['ciii',string(ciii.wvi)],$
            ['caii',string(caii.wvi)],$
            ['mgii',string(mgii.wvi)],$
            ['fkii',string(fkii.wvi)],$
            ['feii',string(feii.wvi)] $
            ]
     info = transpose(info)
     wrest = float(info[*,1])
     mn = min(dblt_name-wrest,imn,/abs)
     return, dblt_retrieve(info[imn,0])
  endelse 

end
