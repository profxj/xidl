;+ 
; NAME:
; ximgd_getct
;    Version 1.0
;
; PURPOSE:
; Read in a pre-defined color table, and invert if necessary.
;
; CALLING SEQUENCE:
;   
;   ximgd_getct, state, flg
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
;   ximgd_getct, state, flg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro ximgd_getct, state, tablenum, CLR=clr

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'ximgd_getct, state, tablenum'
    return
  endif 

  common xcommon_color


  loadct, tablenum, /silent, bottom=8
  if keyword_set(CLR) then tvlct, r, g, b, 8, /get else $
    tvlct, r, g, b, /get

  ; Init colors
  ximgd_initcolors, CLR=clr

  r = r[0:state.ncolors-2]
  g = g[0:state.ncolors-2]
  b = b[0:state.ncolors-2]

  r_vector = r
  g_vector = g
  b_vector = b

  ximgd_stretchct, state

  return
end
