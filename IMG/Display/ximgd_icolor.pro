;+ 
; NAME:
; ximgd_icolor
;    Version 1.1
;
; PURPOSE:
; Sets the color for the Image given a name or index
;
; CALLING SEQUENCE:
;   ximgd_icolor, color
;
; INPUTS:
;   color  -- Color index (string or long)
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
;   ximgd_icolor, 'red'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP (taken from some Schlegel program)
;-
;------------------------------------------------------------------------------


pro ximgd_icolor, color

if (n_elements(color) EQ 0) then return, 1

ncolor = N_elements(color)

; If COLOR is a string or array of strings, then convert color names
; to integer values
if (size(color,/tname) EQ 'STRING') then begin ; Test if COLOR is a string
    
; Detemine the default color for the current device
    if (!d.name EQ 'X') then defcolor = 7 $ ; white for X-windows
    else defcolor = 0           ; black otherwise
    
    icolor = 0 * (color EQ 'black') $
      + 1 * (color EQ 'red') $
      + 2 * (color EQ 'green') $
      + 3 * (color EQ 'blue') $
      + 4 * (color EQ 'cyan') $
      + 5 * (color EQ 'magenta') $
      + 6 * (color EQ 'yellow') $
      + 7 * (color EQ 'white') $
      + defcolor * (color EQ 'default')
    
endif else begin
    icolor = long(color)
endelse

return, icolor
end 
