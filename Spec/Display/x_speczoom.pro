;+ 
; NAME:
;  x_speczoom   Version 1.0
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

pro x_speczoom, state, flg

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_speczoom, state, flg [V1.0]'
    return
  endif 

;  Optional Keywords

  case flg of
      0: begin ; zoom in 
          ; Set center
          center = xgetx_plt(state, /strct)
          deltx = (state.xymnx[2]-state.xymnx[0])/4.

          state.xymnx[0] = center - deltx
          state.xymnx[2] = center + deltx
      end
      1: begin ; zoom out
          center = xgetx_plt(state, /strct)
          deltx = state.xymnx[2]-state.xymnx[0]

;          state.xymnx[0] = (center - deltx) > state.svxymnx[0]
;          state.xymnx[2] = (center + deltx) < state.svxymnx[2]
          state.xymnx[0] = (center - deltx) 
          state.xymnx[2] = (center + deltx) 
      end
      else :
  endcase
  return
end
