;+ 
; NAME:
; xdimg_delov
;     Version 1.1
;
; PURPOSE:
;   Delete a set of OV images
;
; CALLING SEQUENCE:
;  xdimg_delov, strct, delimg
;
; INPUTS:
;   strct  -- Direct image structure
;   delimg -- Indices of direct images
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
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro xdimg_delov, struct, delimg

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'xdimg_delov, struct, delimg [v1.1]'
      return
  endif 

;
  ndel = n_elements(delimg)
  print, 'Deleting the ov files'
  for q=0,ndel-1 do spawn, '\rm '+struct[delimg[q]].img_ov
  struct[delimg].flg_ov = 0

  return
end
  
      
