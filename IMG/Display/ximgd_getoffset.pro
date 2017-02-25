;+ 
; NAME:
; ximgd_getoffset
;    Version 1.1
;
; PURPOSE:
; Routine to calculate the display offset for the current value of
; state.centerpix, which is the central pixel in the display window.
;
; CALLING SEQUENCE:
;   
;   ximgd_getoffset, state
;
; INPUTS:
;   state       - Structure with tv dependent info
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
;   ximgd_getoffset, state
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro ximgd_getoffset, state

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'ximgd_getoffset, state'
    return
  endif 


  state.zoom.offset = round( state.zoom.centerpix - $
                             (0.5 * state.tv.winsize / state.zoom.factor) )
  return
end
