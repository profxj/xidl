pro fill_elmxh, stddla

  close, /all
  dumc = ''
  dumi = 0

; Create the abundance array

  ndla = n_elements(stddla)

  for j=0,ndla-1 do begin

      openr, 2, slc(stddla[j].abndfil)
      readf, 2, dumc
      readf, 2, dumi
      stddla[j].ndfil = dumi
      if(dumi mod 2 GT 0) then begin
          readf, 2, dumc 
          stddla[j].Hfil = dumc
      endif
      if(dumi mod 4 GT 1) then begin
          readf, 2, dumc 
          stddla[j].Efil = dumc
      endif
      if(dumi mod 8 GT 3) then begin
          readf, 2, dumc 
          stddla[j].Ufil = dumc
      endif
      if(dumi mod 16 GT 7) then begin
          readf, 2, dumc 
          stddla[j].Xfil = dumc
      endif
      readf, 2, dumc
      readf, 2, dumc
      XHfil = strmid(dumc,0,xlc(dumc)-3)+'XH'
      
      readfmt, XHfil, 'i2,1x,f6.3,2x,f5.3,1x,i1,2x,i2',$
        dumi1, dumr1, dumr2, dumi2, dumi3, /silent
      for k=0,n_elements(dumi1)-1 do begin
; elm
          stddla[j].elm[dumi1[k]].clm = dumr1[k] 
          stddla[j].elm[dumi1[k]].sigclm = dumr2[k]
          stddla[j].elm[dumi1[k]].flgclm = dumi2[k]
; XH
          getabnd, nm, dumi1[k], abnd, flag=1
          stddla[j].XH[dumi1[k]].clm = dumr1[k] - stddla[j].NHI + 12.0 - abnd
          if stddla[j].sigNHI[1] GT 0. then sigNHI = mean(stddla[j].sigNHI) else $
            sigNHI = stddla[j].sigNHI[0]
          stddla[j].XH[dumi1[k]].sigclm = sqrt(dumr2[k]^2 + sigNHI^2)
          stddla[j].XH[dumi1[k]].flgclm = dumi2[k]
      endfor
      close, 2
  endfor

return
end
