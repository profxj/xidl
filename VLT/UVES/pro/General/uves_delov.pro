;+ 
; NAME:
; uves_delov   
;     Version 1.0
;
; PURPOSE:
;    Deletes OV images
;
; CALLING SEQUENCE:
;   
;  uves_delov, uves, indx
;
; INPUTS:
;   uves   -  ESI structure
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
;   uves_delov, uves, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro uves_delov, uves, indx, SILENT=silent

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'uves_delov, uves, indx, /silent (v1.0)'
      return
  endif 

;
  ndel = n_elements(indx)
  if not keyword_set(SILENT) then print, 'uves_delov: Deleting the ov files'
  for q=0L,ndel-1 do begin
      outfil = uves_getfil('ov_fil', $
                            OBJN=uves[indx[q]].img_root, /name)
      if not keyword_set(SILENT) then print, '\rm '+outfil
      spawn, '\rm '+outfil
  endfor

  uves[indx].flg_ov = 0

  return
end
  
      
