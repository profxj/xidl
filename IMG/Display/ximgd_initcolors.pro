;+ 
; NAME:
; ximgd_initcolors
;    Version 1.0
;
; PURPOSE:
; Sets up zooming
;
; CALLING SEQUENCE:
;   
;   ximgd_initcolors, state, flg
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
;   ximgd_initcolors, state, flg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro ximgd_initcolors, CLR=clr

; Load a simple color table with the basic 8 colors in the lowest 
; 8 entries of the color table.  Also set top color to white.

  if keyword_set( CLR ) then begin
      rtiny   = [0, 1, 0, 0, 0, 1, 1, 1]
      gtiny = [0, 0, 1, 0, 1, 0, 1, 1]
      btiny  = [0, 0, 0, 1, 1, 1, 0, 1]
      tvlct, 255*rtiny, 255*gtiny, 255*btiny
  endif

  tvlct, [255],[255],[255], !d.table_size-1

end
