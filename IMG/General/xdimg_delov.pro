pro xdimg_delov, struct, delimg

  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'xdimg_delov, struct, delimg (v1.0)'
      return
  endif 

;
  ndel = n_elements(delimg)
  print, 'Deleting the ov files'
  for q=0,ndel-1 do spawn, '\rm '+struct[delimg[q]].img_ov
  struct[delimg].flg_ov = 0

  return
end
  
      
