;+ 
; NAME:
; fill_ion
;  V1.0
;
; PURPOSE:
;    Fills up ions
;
; CALLING SEQUENCE:
;   
;   fill_ion, sDLA
;
; INPUTS:
;
; RETURNS:
;   structure      - IDL structure
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
;   fill_ion, sdla
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   12-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------
pro fill_ion, stddla

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'fill_ion, struct (v1.0)' 
    return
  endif 
;

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
      readf, 2, dumc
      readf, 2, dumc
      close, 2
; Ionfil
      ionfil = strmid(dumc,0,xlc(dumc)-3)+'all'
      if numlines(ionfil) EQ 0 then continue

      readcol, ionfil, format='I,I,F,F,I,I', dumi1, dumi2, dumr1, $
        dumr2, dumi3, dumi4, /silent

; ion
      for k=0,n_elements(dumr1)-1 do begin
          stddla[j].ion[dumi1[k]].state[dumi2[k]].clm = dumr1[k]
          stddla[j].ion[dumi1[k]].state[dumi2[k]].sigclm = dumr2[k]
          stddla[j].ion[dumi1[k]].state[dumi2[k]].flgclm = dumi3[k]
      endfor

  endfor

return
end
