;+ 
; NAME:
; hires_delov   
;     Version 1.0
;
; PURPOSE:
;    Deletes OV images
;
; CALLING SEQUENCE:
;  hires_delov, hires, indx
;
; INPUTS:
;   hires   -  ESI structure
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
;   hires_delov, hires, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro hires_delov, hires, indx, SILENT=silent

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_delov, hires, indx, /silent (v1.0)'
      return
  endif 

;
  ndel = n_elements(indx)
  if not keyword_set(SILENT) then print, 'hires_delov: Deleting the ov files'
  for q=0L,ndel-1 do begin
      outfil = hires_getfil('ov_fil', $
                            CHIP=hires[indx[q]].chip, $
                            FRAME=hires[indx[q]].frame, /name)
      if not keyword_set(SILENT) then print, '\rm '+outfil
      spawn, '\rm '+outfil
  endfor

  hires[indx].flg_ov = 0

  return
end
  
      
