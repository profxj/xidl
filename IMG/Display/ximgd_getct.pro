;+ 
; NAME:
; ximgd_getct
;    Version 1.1
;
; PURPOSE:
; Read in a pre-defined color table
;
; CALLING SEQUENCE:
;   
;   ximgd_getct, state, tablenum, /CLR
;
; INPUTS:
;   state       - Structure with tv dependent info
;   tablenum    - Index of color table
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /CLR -- Sets to 8 bit?
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   ximgd_getct, state, 1
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro ximgd_getct, state, tablenum, CLR=clr, NOTRIM=notrim

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

  if not keyword_set(NOTRIM) then begin
     r = r[0:state.ncolors-2]
     g = g[0:state.ncolors-2]
     b = b[0:state.ncolors-2]
  endif

  r_vector = r
  g_vector = g
  b_vector = b

  ximgd_stretchct, state

  return
end
