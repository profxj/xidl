;+ 
; NAME:
; esi_delov   
;     Version 1.0
;
; PURPOSE:
;    Deletes OV images
;
; CALLING SEQUENCE:
;   
;  esi_delov, esi, indx
;
; INPUTS:
;   esi   -  ESI structure
;   indx  -  Index numbers of frames to delete
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
;   esi_delov, esi, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro esi_delov, esi, indx, SILENT=silent

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_delov, esi, indx, /silent (v1.0)'
      return
  endif 

;
  ndel = n_elements(indx)
  if not keyword_set(SILENT) then print, 'esi_delov: Deleting the ov files'
  for q=0L,ndel-1 do spawn, '\rm '+esi[indx[q]].img_ov
  esi[indx].flg_ov = 0

  return
end
  
      
