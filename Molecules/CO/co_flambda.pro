;+ 
; NAME:
; co_flambda
;  (V1.0)
;
; PURPOSE:
;    Calculate the flambda value for a CO A-X transition
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

;; This assumes an N-0  A-X transition (i.e. vpp=0)
function co_flambda, vp, jp, jpp  


  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'flambda = co_flambda(vp, jp, jpp) [v1.0]'
    return, -1
  endif 
  ;; log flambda
  ;; Normalize to R(0) levels of Morton
  ;; There may be a better choice (less perturbed)
  case vp of
      0: nrm_fl = 1.387
      1: nrm_fl = 1.643
      2: nrm_fl = 1.773
      3: nrm_fl = 1.700
      4: nrm_fl = 1.520
      5: nrm_fl = 1.305
      6: nrm_fl = 1.036
      7: nrm_fl = 0.745
      8: nrm_fl = 0.423
      9: nrm_fl = 0.092
      else: stop
  endcase

  ;; H-L factor for R(0) is unity, i.e. ignore

  ;; H-L term
  case jp-jpp of 
      -1: HLfact = float(jpp-1)/(2*(2*jpp+1)) ;; P
      0: HLfact = 0.5 ;; Q
      1: HLfact = float(jpp+2)/(2*(2*jpp+1)) ;; R
      else: stop
  endcase

  ;; Scale
  flamb =  10.^nrm_fl * HLfact  ; Ang

  ;;
  return, flamb
end
