;+ 
; NAME:
; apf_delov   
;     Version 1.0
;
; PURPOSE:
;    Deletes OV images
;
; CALLING SEQUENCE:
;  apf_delov, apf, indx
;
; INPUTS:
;   apf   -  ESI structure
;   indx  -  Index numbers of frames to delete
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; /SILENT -- No messages printed to screen
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   apf_delov, apf, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro apf_delov, apf, indx, SILENT=silent

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'apf_delov, apf, indx, /silent (v1.0)'
      return
  endif 

;
  ndel = n_elements(indx)
  if not keyword_set(SILENT) then print, 'apf_delov: Deleting the ov files'
  for q=0L,ndel-1 do begin
      outfil = apf_getfil('ov_fil',  FRAME=apf[indx[q]].frame, /name)
      if not keyword_set(SILENT) then print, '\rm '+outfil
      spawn, '\rm '+outfil
  endfor

  apf[indx].flg_ov = 0

  return
end
  
      
