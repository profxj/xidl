;+ 
; NAME:
; h2_xlevel
;  (V1.0)
;
; PURPOSE:
;    Returns the energy relative to true ground for H2 in the
;    electronic ground-state X for a given J value
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;  energy (cm^-1) relative to the true ground state (v=0, J=0)
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
function h2_xlevel, jpp

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'engy = h2_xlevel(jpp) [v1.0]'
    return, -1
  endif 
  ;; Values were given by Sheffer
  case jpp of
      0: engy=    0.000d
      1: engy=  118.483d
      2: engy=  354.366d
      3: engy=  705.511d
      4: engy= 1168.791d
      5: engy= 1740.178d
      6: engy= 2414.856d
      7: engy= 3187.355d
      8: engy= 4051.685d
      9: engy= 5001.488d
      10: engy= 6030.171
      11: engy= 7131.019d
      12: engy= 8297.289d
      13: engy= 9522.243d
      else: stop
  endcase
  ;;
  return, engy  ; cm^-1
end
