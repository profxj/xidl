;+ 
; NAME:
; esi_delfin   
;     Version 1.0
;
; PURPOSE:
;    Deletes Final images
;
; CALLING SEQUENCE:
;   
;  esi_delfin, esi, indx
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
;   esi_delfin, esi, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro esi_delfin, esi, indx, SILENT=silent

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_delfin, esi, indx, /silent (v1.0)'
      return
  endif 

;
  ndel = n_elements(indx)
  if not keyword_set(SILENT) then print, 'esi_delfin: Deleting the ov files'
  for q=0L,ndel-1 do spawn, '\rm '+esi[indx[q]].img_final+'.gz'
  esi[indx].flg_ov = 0

  return
end
  
      
