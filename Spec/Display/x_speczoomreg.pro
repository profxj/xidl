;+ 
; NAME:
;  x_speczoomreg
;    Version 1.1
;
; PURPOSE:
;    Sets region to zoom in on in a GUI
;
; CALLING SEQUENCE:
;   x_speczoom, state
;
; INPUTS:
;  state  -- GUI state (requires TAGS xymnx tmpxy, flg_zoom)
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
;   x_speczoom, state
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_speczoomreg, state

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_speczoomreg, state [v1.1]'
    return
  endif 

;  Optional Keywords

  if state.flg_zoom EQ 0 then begin
          state.tmpxy[0] = xgetx_plt(state, /strct)
          state.flg_zoom = 1
  end else begin
      state.tmpxy[2] = xgetx_plt(state, /strct)
      state.flg_zoom = 2
      tmpxy = state.xymnx
      state.xymnx[0] = state.tmpxy[0] < state.tmpxy[2]
      state.xymnx[2] = state.tmpxy[0] > state.tmpxy[2]
  endelse

  return
end
