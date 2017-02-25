;+ 
; NAME:
; hires_slitlen
;     Version 1.1
;
; PURPOSE:
;  Returns the slitlenght (arcseconds) given the decker name
;
; CALLING SEQUENCE:
;  slitlen = hires_slitlen(decker)
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
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Aug-2005 Written by JXP  
;-
;------------------------------------------------------------------------------

function hires_slitlen, decker

  ;;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'slit_len = hires_slitlen(decker) [v1.0]'
      return, -1
  endif 

  case strtrim(decker,2) of
      'B1': slit_len = 3.5
      'B2': slit_len = 7.0
      'B3': slit_len = 14.0
      'B4': slit_len = 28.0
      'B5': slit_len = 3.5
      'C1': slit_len = 7.0
      'C2': slit_len = 14.0
      'C3': slit_len = 28.0
      'C4': slit_len = 3.5
      'C5': slit_len = 7.0
      'D1': slit_len = 14.0
      'D2': slit_len = 28.0
      'D3': slit_len = 7.
      'D4': slit_len = 14.
      'E1': slit_len = 5.
      'E3': slit_len = 5.
      'E4': slit_len = 7.
      else: stop  ;; input slit_len by hand and inform JXP 
  endcase

  return, slit_len
end

