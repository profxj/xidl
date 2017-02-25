;+ 
; NAME:
;  x_specpan   
;   Version 1.1
;
; PURPOSE:
;    Moves about a spectrum to left or right in a GUI.  The
;  state structure must have TAGS: xymnx and wave
;
; CALLING SEQUENCE:
;   x_specpan, state, /LEFT, /NOY
;
; INPUTS:
;  state -- GUI state structure with TAGS: xymnx, wave
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /LEFT -- Pan left (default is right)
;  /NOY  -- Do not rescale in y  (default is to rescale)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_specpan, state, /LEFT
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_specpan, state, LEFT=left, NOY=noy
;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_specpan, state, /LEFT, /NOY [V1.1]'
    return
  endif 

  ; Size of screen
  sz = state.xymnx[2]-state.xymnx[0]

  ; Right or left
  if not keyword_set( LEFT ) then begin ; RIGHT
      state.xymnx[0] = state.xymnx[2] - 0.1*sz
      state.xymnx[2] = state.xymnx[2] + 0.9*sz
  endif else begin
      state.xymnx[2] = state.xymnx[0] + 0.1*sz
      state.xymnx[0] = state.xymnx[0] - 0.9*sz
  endelse

  ; Reset top and bottom
  if not keyword_set( NOY ) then begin
      gdwv = where(state.wave GT state.xymnx[0] AND $
                   state.wave LT state.xymnx[2], ngd)
      if ngd GT 0 then begin
          mn = min(state.fx[gdwv], max=mx)
          state.xymnx[1] = mn - (mx-mn)*0.05
          state.xymnx[3] = mx + (mx-mn)*0.05
      endif
  endif

  return
end
