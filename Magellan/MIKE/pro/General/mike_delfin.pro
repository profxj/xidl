;+ 
; NAME:
; mike_delfin   
;     Version 1.0
;
; PURPOSE:
;    Deletes Final images
;
; CALLING SEQUENCE:
;   
;  mike_delfin, mike, indx
;
; INPUTS:
;   mike   -  ESI structure
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
;   mike_delfin, mike, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro mike_delfin, mike, indx, SILENT=silent

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_delfin, mike, indx, /silent (v1.0)'
      return
  endif 

;
  ndel = n_elements(indx)
  if not keyword_set(SILENT) then print, 'mike_delfin: Deleting the ov files'
  for q=0L,ndel-1 do spawn, '\rm '+mike[indx[q]].img_final+'.gz'
  mike[indx].flg_final = 0

  return
end
  
      
