;+ 
; NAME:
; ximgd_contrast
;    Version 1.1
;
; PURPOSE:
;   Routine to set the contrast and brightness depending on the 
;   position of the cursor.
;
; CALLING SEQUENCE:
;   ximgd_contrast, state
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
;   ximgd_contrast, state
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro ximgd_contrast, state

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'ximgd_contrast, state'
    return
  endif 

  common xcommon_color

  state.contrast = state.tv.xcurs/float(state.tv.winsize[0])
  state.brightness = state.tv.ycurs/float(state.tv.winsize[1])
  ximgd_stretchct, state
  state.zoom.flg = 2

  return

end
