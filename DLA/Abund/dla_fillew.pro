;+ 
; NAME:
; dla_fillew
;  V1.0
;
; PURPOSE:
;  Fills up the arrays with EW values using the .EW files
;
; CALLING SEQUENCE:
;
; INPUTS:
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
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2006 Written by JXP
;-
;------------------------------------------------------------------------------
pro dla_fillew, stddla, ROOT=root

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'dla_fillew, dla [v1.0]' 
    return
  endif 

  if not keyword_set(ROOT) then root = ''
;

  close, 2
  dumc = ''
  dumi = 0

; Create the abundance array

  ndla = n_elements(stddla)

  for j=0,ndla-1 do begin

      ;; Abund fil
      a = file_search(root+strtrim(stddla[j].abndfil,2),count=na)
      if na EQ 0 then begin
          print, 'dla_fillew: No Abund file -- ', root+stddla[j].abndfil
          continue
      endif

      openr, 2, strtrim(root+stddla[j].abndfil,2)
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
      if(dumi mod 128 GT 63) then begin
          readf, 2, dumc 
          stddla[j].Ffil = dumc
      endif
      readf, 2, dumc
      readf, 2, dumc
      close, 2


      ;; Allfil
      ionfil = root+strmid(dumc,0,strlen(strtrim(dumc,2))-3)+'all'
      a = file_search(ionfil,count=na)
      if na EQ 0 then begin
         print, 'dla_fillew: No All file -- ', ionfil
         continue
      endif else begin
         if numlines(ionfil) EQ 0 then begin 
            print, 'dla_fillew: Empty file-- ', ionfil
            continue
         endif
      endelse

      readcol, ionfil, format='I,I,F,F,I,I', dumi1, dumi2, dumr1, $
        dumr2, dumi3, dumi4, /silent

      ;; Central ion
      for k=0,n_elements(dumr1)-1 do begin
          stddla[j].ion[dumi1[k]].state[dumi2[k]].clm = dumr1[k]
          stddla[j].ion[dumi1[k]].state[dumi2[k]].sigclm = dumr2[k]
          stddla[j].ion[dumi1[k]].state[dumi2[k]].flgclm = dumi3[k]
      endfor

      ;; EW fil
      ionfil = root+strmid(dumc,0,strlen(strtrim(dumc,2))-3)+'EW'
      a = findfile(ionfil,count=na)
      if na EQ 0 then begin
          print, 'dla_fillew: No EW file -- ', ionfil
          continue
      endif
      if numlines(ionfil) EQ 0 then continue

      readcol, ionfil, format='D,F,F,I,I', wrest, EW, sigEW, inst, /sile

      ;; Central ion
      for k=0,n_elements(wrest)-1 do begin
          ;; Grab atomic number and ion
          getion, wrest[k], ion, Z=z
          stddla[j].ion[Z].indx[ion]++
          idx = stddla[j].ion[Z].indx[ion]
          stddla[j].ion[Z].state[ion,idx].clm = EW[k]      ; Rest EW
          stddla[j].ion[Z].state[ion,idx].lambda = wrest[k]
          stddla[j].ion[Z].state[ion,idx].sigclm = sigEW[k]
;          stddla[j].ion[Z].state[ion,idx].flgclm = dumi1[k]
          stddla[j].ion[Z].state[ion,idx].flginst = inst[k]
      endfor

  endfor

return
end
