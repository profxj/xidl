;+ 
; NAME:
;  x_speczoom   
;   Version 1.1
;
; PURPOSE:
;    Zooms in or out on a spectrum in a GUI
;
; CALLING SEQUENCE:
;   x_speczoom, state, flg
;
; INPUTS:
;   state - GUI state structure (needs TAG xymnx)
;   flg   - 0-In, 1-Out
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
;  xgetx_plt
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
             'x_speczoom, state, flg [v1.1]'
    return
  endif 

;  Optional Keywords

  if tag_exist(state, 'flg_logx', /top_level) then log = state.flg_logx $
    else log = 0

  case flg of
      0: begin ; zoom in 
          ; Set center
          center = xgetx_plt(state, /strct)
          if LOG then begin
              deltx = (alog10(state.xymnx[2])-alog10(state.xymnx[0]))/4.
              state.xymnx[0] = 10.^(alog10(center) - deltx)
              state.xymnx[2] = 10.^(alog10(center) + deltx)
          endif else begin
              deltx = (state.xymnx[2]-state.xymnx[0])/4.
              state.xymnx[0] = center - deltx
              state.xymnx[2] = center + deltx
          endelse
      end
      1: begin ; zoom out
          center = xgetx_plt(state, /strct)
          if LOG then begin
              deltx = alog10(state.xymnx[2])-alog10(state.xymnx[0])
              state.xymnx[0] = 10.^(alog10(center) - deltx) 
              state.xymnx[2] = 10.^(alog10(center) + deltx) 
          endif else begin
              deltx = state.xymnx[2]-state.xymnx[0]
              state.xymnx[0] = (center - deltx) 
              state.xymnx[2] = (center + deltx) 
          endelse
      end
      else :
  endcase
  return
end
