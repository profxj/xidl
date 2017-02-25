;+ 
; NAME:
; ximgd_stretchct
;    Version 1.1
;
; PURPOSE:
; routine to change color stretch for given values of 
; brightness and contrast.
; This routine is now shorter and easier to understand.  
; It is based on Baarth and Finkbeinder code.
;
; CALLING SEQUENCE:
;   
;   ximgd_stretchct, state
;
; INPUTS:
;   state       - Structure with image info (e.g. brightness, contrast)
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
;   ximgd_stretchct, state
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
; Complete rewrite 2000-Sep-21 - Doug Finkbeiner
;   08-Feb-2002 Revised by JXP
;-
;------------------------------------------------------------------------------


pro ximgd_stretchct, state

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'ximgd_stretchct, state'
    return
  endif 

  common xcommon_color

  x = state.brightness*(state.ncolors-1)
  y = state.contrast*(state.ncolors-1) > 2 ; Minor change by AJB 
  high = x+y & low = x-y
  diff = (high-low) > 1
  
  slope = float(state.ncolors-1)/diff ;Scale to range of 0 : nc-1
  intercept = -slope*low
  p = long(findgen(state.ncolors)*slope+intercept) ;subscripts to select
  tvlct, r_vector[p], g_vector[p], b_vector[p], 8

end

