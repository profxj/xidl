;+ 
; NAME:
; xicmm_colors    
;  Version 1.1
;
; PURPOSE:
;   Initializes the common block for Image colors.
;   For some reason the colors are deleted.
;
; CALLING SEQUENCE:
;   xicmm_colors
;
; INPUTS:
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
;   xicmm_colors
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro xicmm_colors

common xcommon_color, r_vector, g_vector, b_vector

;  delvarx, r_vector, g_vector, b_vector
  r_vector=0 
  g_vector=0 
  b_vector=0

end
