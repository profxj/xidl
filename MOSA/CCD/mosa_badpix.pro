pro mosa_badpix, finlst, CLFIL=clfil, IM3=im3, IBAND=iband

  if not keyword_set( clfil ) then clfil = 'badpix.cl'

  ;; Read list
  readcol, finlst, fil, format='A'

  close, /all
  openw, 1, clfil
  ;; Loop
  for q=0L,n_elements(fil)-1 do begin
      ;; IM1
      if keyword_set( IBAND ) then $
      printf, 1, 'imreplace bpm'+strtrim(strmid(fil[q],3),2)+ $
        '/bpm_im1.pl[82:163,2020:2108] 1'
      ;; IM3
      if keyword_set( IM3) then $
        printf, 1, 'imreplace bpm'+strtrim(strmid(fil[q],3),2)+ $
        '/bpm_im3.pl[1353:1354,1:1986] 1'
  endfor

  close, /all
  return
end

  
