pro parse_DLAion, strct

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'parse_DLAion, dlastruct'
    return
  endif 

;  Loop on the strcuture

  for i=0,n_elements(struct) do begin
      
;  Open the Ion file

      

