;+ 
; NAME:
;  x_speczoomreg
;    Version 1.0
;
; PURPOSE:
;    Sets a wavelength array given a header
;
; CALLING SEQUENCE:
;   
;   x_speczoom, state, flg
;
; INPUTS:
;   flg  = 0-In, 1-Out
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
;   x_speczoom, state, flg
;
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
             'x_speczoomreg, state [V1.0]'
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
