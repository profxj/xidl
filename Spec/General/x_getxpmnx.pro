;+ 
; NAME:
; x_getxpmnx
;
; PURPOSE:
;  Create an array of velocity arrays for a string of transitions
;
; CALLING SEQUENCE:
;   
;   all_velo = x_getxpmnx(wave, zabs, wrest, vmnx, ALL_PMNX=)
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
;   all_velo = x_getxpmnx(wave, zabs, wrest, vmnx, ALL_PMNX=)
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   24-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------

function x_getxpmnx, state

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'xpmnx = x_getxpmnx(state) [v1.0]'
    return, -1
  endif 

  a = min(abs(state.wave - state.xymnx[0]), imn)
  b = min(abs(state.wave - state.xymnx[2]), imx)

  return, [imn,imx]

end
  
