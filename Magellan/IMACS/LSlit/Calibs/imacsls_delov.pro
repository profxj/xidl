;+ 
; NAME:
; imacsls_delov   
;     Version 1.1
;
; PURPOSE:
;    Deletes OV images
;
; CALLING SEQUENCE:
;  imacsls_delov, imacsls, indx
;
; INPUTS:
;   imacsls -  IMACS structure
;   indx    -  Index numbers of frames to delete
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
;   imacsls_delov, imacsls, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro imacsls_delov, imacsls, indx, SILENT=silent

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'imacsls_delov, imacsls, indx, /silent [v1.1]'
      return
  endif 

;
  ndel = n_elements(indx)
  if not keyword_set(SILENT) then print, 'imacsls_delov: Deleting the ov files'
  for q=0L,ndel-1 do spawn, '\rm '+imacsls[indx[q]].img_ov
  imacsls[indx].flg_ov = 0

  return
end
  
      
