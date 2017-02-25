;+ 
; NAME:
; ximgd_initcolors
;    Version 1.1
;
; PURPOSE:
; Load a simple color table with the basic 8 colors in the lowest 
; 8 entries of the color table.  Also set top color to white.
;
; CALLING SEQUENCE:
;   
;   ximgd_initcolors, /CLR
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /CLR  -- Sets up the colors (not sure why you wouldnt set CLR)
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
