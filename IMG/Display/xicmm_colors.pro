;+ 
; NAME:
; xicmm_colors    
;  Version 1.0
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
;   xicmm_colors, state, flg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro xicmm_colors

common xcommon_color, r_vector, g_vector, b_vector

  delvarx, r_vector, g_vector, b_vector

end
