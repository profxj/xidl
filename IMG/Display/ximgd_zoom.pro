;+ 
; NAME:
; ximgd_zoom
;    Version 1.1
;
; PURPOSE:
;  Change the zoom level and recenter as desired
;
; CALLING SEQUENCE:
;   ximgd_zoom, state, inout, /recenter
;
; INPUTS:
;   state       - Structure with tv dependent info
;   inout       - string ('in' or 'out')
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  /RECENTER  -- Recenter on cursor position.
;
; COMMENTS:
;
; EXAMPLES:
;   ximgd_zoom, state, 'in', /recenter
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro ximgd_zoom, state, inout, RECENTER=recenter

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'ximgd_zoom, state, inout, /recenter'
    return
  endif 

  case inout of
      'in':    state.zoom.level = (state.zoom.level + 1) < 12
      'out':   state.zoom.level = (state.zoom.level - 1) > (-12) 
      'one':   state.zoom.level =  0
      'none':  ; no change to zoom level: recenter on current mouse position
      else:  print,  'problem in xatv_zoom!'
  endcase

  ; ZOOMLEVEL
  if state.zoom.level LE 0 then $
    state.zoom.factor = 1./(1.+abs(state.zoom.level)*0.5) else $
    state.zoom.factor = 2.^state.zoom.level
  
  ; RECENTER
  if keyword_set( RECENTER ) then begin
      state.zoom.centerpix[0] = x_tvx(state.tv, /intg)
      state.zoom.centerpix[1] = x_tvy(state.tv, /intg)
      ximgd_getoffset, state
  endif


  if keyword_set( RECENTER ) then begin
      newpos = (state.zoom.centerpix - state.zoom.offset + 0.5) * state.zoom.factor
      tvcrs, newpos[0], newpos[1], /device 
  endif

  return

end
