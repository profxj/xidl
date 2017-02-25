;+ 
; NAME:
; x_getxpmnx
;
; PURPOSE:
;  Find the pixels corresponding to xymnx in state
;
; CALLING SEQUENCE:
;   pxmnx = x_getxpmnx(state)
;
; INPUTS:
;  state -- Structure with TAGS wave, xymnx
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
;   pxmnx = x_getxpmnx(state)
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
             'xpmnx = x_getxpmnx(state) [v1.1]'
    return, -1
  endif 

  a = min(abs(state.wave - state.xymnx[0]), imn)
  b = min(abs(state.wave - state.xymnx[2]), imx)

  return, [imn,imx]

end
  
