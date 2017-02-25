;+ 
; NAME:
; hamspec_delov   
;     Version 1.0
;
; PURPOSE:
;    Deletes OV images
;
; CALLING SEQUENCE:
;  hamspec_delov, hamspec, indx
;
; INPUTS:
;   hamspec   -  ESI structure
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
;   hamspec_delov, hamspec, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro hamspec_delov, hamspec, indx, SILENT=silent

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_delov, hamspec, indx, /silent (v1.0)'
      return
  endif 

;
  ndel = n_elements(indx)
  if not keyword_set(SILENT) then print, 'hamspec_delov: Deleting the ov files'
  for q=0L,ndel-1 do begin
      outfil = hamspec_getfil('ov_fil', $
                            FRAME=hamspec[indx[q]].frame, /name)
      if not keyword_set(SILENT) then print, '\rm '+outfil
      spawn, '\rm '+outfil
  endfor

  hamspec[indx].flg_ov = 0

  return
end
  
      
