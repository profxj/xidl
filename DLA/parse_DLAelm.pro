pro parse_DLAelm, supstrc, stddla

  close, /all
  dumc = ''
  dumi = 0
  ndla = n_elements(stddla)

  strctnm = {elmDLA, dla_flgelm: intarr(100),$
	dla_elm: fltarr(100),$
	dla_elmsig: fltarr(100) }
	
; Create the abundance array
  supstrc = replicate(strctnm,ndla)

  for j=0,ndla-1 do begin
    	openr, 2, slc(stddla[j].dla_abndfil)
	readf, 2, dumc
	readf, 2, dumi
	stddla[j].dla_ndfil = dumi
   	if(dumi mod 2 GT 0) then begin
	  readf, 2, dumc 
	  stddla[j].dla_Hfil = dumc
	endif
   	if(dumi mod 4 GT 1) then begin
	  readf, 2, dumc 
	  stddla[j].dla_Efil = dumc
	endif
   	if(dumi mod 8 GT 3) then begin
	  readf, 2, dumc 
	  stddla[j].dla_Ufil = dumc
	endif
	readf, 2, dumc
	readf, 2, dumc
	XHfil = strmid(dumc,0,xlc(dumc)-3)+'XH'

	readfmt, XHfil, 'i2,1x,f6.3,2x,f5.3,1x,i1,2x,i2',$
		dumi1, dumr1, dumr2, dumi2, dumi3, /silent
	for k=0,n_elements(dumi1)-1 do begin
	  supstrc[j].dla_elm[dumi1[k]] = dumr1[k]
	  supstrc[j].dla_elmsig[dumi1[k]] = dumr2[k]
	  supstrc[j].dla_flgelm[dumi1[k]] = dumi2[k]
	endfor
  	close, 2
  endfor

return
end
