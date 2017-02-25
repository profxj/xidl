;+ 
; NAME:
; xdimg_deldrksub
;     Version 1.1
;
; PURPOSE:
;   Delete a set of dark subtracted images.
;
; CALLING SEQUENCE:
;  xdimg_deldrksub, strct, delimg
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
;   14-June-2007 Written by LKP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro xdimg_deldrksub, struct, delimg

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'xdimg_deldrksub, struct, delimg [v1.1]'
      return
  endif 

  ndel = n_elements(delimg)
  print, 'Deleting the dark subtracted files.'

  for q=0,ndel-1 do begin
     if struct[delimg[q]].flg_drksub EQ 0 then begin
        print, 'Cannot delete non-existent file.  Skipping!'
        continue
     endif
     spawn, '\rm '+struct[delimg[q]].img_drksub
     struct[delimg[q]].flg_drksub = 0
  endfor
     
  return
end
  
      
