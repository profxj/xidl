pro wfccd_delov, struct, delimg, SILENT=silent

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'wfccd_delov, struct, delimg, /silent (v1.0)'
      return
  endif 

;
  ndel = n_elements(delimg)
  if not keyword_set(SILENT) then print, 'wfccd_delov: Deleting the ov files'
  for q=0L,ndel-1 do spawn, '\rm '+struct[delimg[q]].img_ov
  struct[delimg].flg_ov = 0

  return
end
  
      
