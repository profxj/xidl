pro writecol, name, v1, v2, FORMAT= fmt


; writecol -- Writes a 2 column ascii file

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'writecol, name, v1, v2, [FORMAT= ]'
    return
  endif 

;

  if not keyword_set( FORMAT ) then    flgfmt    = 0 else flgfmt = 1

;

  close, 9
  openw, 9, name
  for i=0,n_elements(v1)-1 do begin
      if (flgfmt = 1) then printf, 9, FORMAT = fmt, v1[i], v2[i] $
      else printf, 9, v1[i], v2[i]
  endfor

  close, 9

return
end
