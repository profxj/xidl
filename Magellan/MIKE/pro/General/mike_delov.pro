;+ 
; NAME:
; mike_delov   
;     Version 1.0
;
; PURPOSE:
;    Deletes OV images
;
; CALLING SEQUENCE:
;   
;  mike_delov, mike, indx
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
;   mike_delov, mike, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro mike_delov, mike, indx, SILENT=silent

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_delov, mike, indx, /silent (v1.0)'
      return
  endif 

;
  ndel = n_elements(indx)
  if not keyword_set(SILENT) then print, 'mike_delov: Deleting the ov files'
  for q=0L,ndel-1 do begin
      outfil =mike_getfil('ov_fil', subfil=mike[indx[q]].img_root, $
                          /name, CHKFIL=chkf)
      print, '\rm '+outfil
      spawn, '\rm '+outfil
  endfor

  mike[indx].flg_ov = 0

  return
end
  
      
