;+ 
; NAME:
; ximgd_resize
;    Version 1.1
;
; PURPOSE:
; Routine to resize the draw window when a top-level resize event
; occurs.
;
; CALLING SEQUENCE:
;   
;   ximgd_resize, state
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
;   ximgd_resize, state
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro ximgd_resize, state

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'ximgd_resize, state'
    return
  endif 


  widget_control, state.base_id, tlb_get_size=tmp_event

  drawpad = (widget_info(state.drawbase_id,/geometry)).xsize - $
    state.tv.winsize[0]

  window = (400L > tmp_event)

  newbase = window 

  widget_control, state.drawbase_id, xsize = newbase[0], ysize = newbase[1]

  newxsize = (widget_info(state.drawbase_id,/geometry)).xsize - drawpad
  newysize = (widget_info(state.drawbase_id,/geometry)).ysize - drawpad

  widget_control, state.draw_id, xsize = newxsize, ysize = newysize

  state.tv.winsize = [newxsize, newysize]

  return

end
