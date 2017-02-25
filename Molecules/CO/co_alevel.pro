;+ 
; NAME:
; co_alevel
;  (V1.0)
;
; PURPOSE:
;    Calculate the energy levels of the CO molecule A electronic level
;
; CALLING SEQUENCE:
;
; INPUTS:
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
;   Sep-2008 Written by JXP with guidance from Y Sheffer
;-
;------------------------------------------------------------------------------

function co_alevel, vp, jp

;  c = x_constants()
  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'engy = co_alevel(vp, jp) [v1.0]'
    return, -1
  endif 

  ;;  Constants from Tilford & Simmons (1972)
  TE=65075.77
  WE=1518.240
  WEXE=-19.4000
  WEYE=7.6584e-1
  WEZE=-1.4117e-1
  WEAE=1.434e-2
  WEBE=-8.051e-4
  WECE=2.36e-5
  WEDE=-2.90e-7
  BE=1.6115
  ALP=-2.3251e-2
  GAM0=1.5911e-3
  GAM1=-5.7160e-4
  GAM2=8.2417e-5
  GAM3=-5.9413e-6
  GAM4=2.1149e-7
  GAM5=-2.991e-9
  DE=7.29e-6
  BET=-1.05e-7

  ;; Basis
  vv = (vp + 0.5)
  jj = jp*(jp + 1.)

  ;; Add em up
  GV=VV*WE+VV^2*WEXE+VV^3*WEYE+VV^4*WEZE+VV^5*WEAE+VV^6*WEBE+ $
     VV^7*WECE+VV^8*WEDE

  BV=BE+VV*ALP+VV^2*GAM0+VV^3*GAM1+VV^4*GAM2+VV^5*GAM3+VV^6*GAM4+VV^7*GAM5
  
  DV=DE+VV*BET

  ;; Is it really a minus sign for DV?
  FJ=BV*JJ-DV*JJ^2

  ;; Final energy (relative to the vpp=0 XS ground-state)

  engy = TE + FJ + GV  ; cm^-1

  return, engy
end
