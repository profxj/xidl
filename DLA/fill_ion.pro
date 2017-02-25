;+ 
; NAME:
; fill_ion
;  V1.1
;
; PURPOSE:
;    Fills up the column densities for all of the ions observed for a
;    given DLA.  It parses the .ion files produced by dla_updabd.
;    This program is unlikely to be called by anything except parse_dlalst.
;
; CALLING SEQUENCE:
;   fill_ion, sDLA
;
; INPUTS:
;
; RETURNS:
;   sDLA      - IDL DLA structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  ROOT=  Path to DLA tree (e.g. '~/DLA/')
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
pro fill_ion, stddla, ROOT=root

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'fill_ion, dla (v1.1)' 
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
      if(dumi mod 64 GT 31) then begin
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
      if numlines(ionfil) EQ 0 then continue

      readcol, ionfil, format='I,I,F,F,I,I', dumi1, dumi2, dumr1, $
        dumr2, dumi3, dumi4, /silent

      ;; Central ion
      for k=0,n_elements(dumr1)-1 do begin
          stddla[j].ion[dumi1[k]].state[dumi2[k]].clm = dumr1[k]
          stddla[j].ion[dumi1[k]].state[dumi2[k]].sigclm = dumr2[k]
          stddla[j].ion[dumi1[k]].state[dumi2[k]].flgclm = dumi3[k]
      endfor

      ;; Ionfil
      ionfil = root+strmid(dumc,0,strlen(strtrim(dumc,2))-3)+'ion'
      if numlines(ionfil) EQ 0 then continue

      readcol, ionfil, format='D,F,F,I,I', dumd1, dumr1, $
        dumr2, dumi1, dumi2, /silent

      ;; Central ion
      for k=0,n_elements(dumr1)-1 do begin
          ;; Grab atomic number and ion
          getion, dumd1[k], ion, Z=z
          stddla[j].ion[Z].indx[ion]++
          idx = stddla[j].ion[Z].indx[ion]
          stddla[j].ion[Z].state[ion,idx].clm = dumr1[k]
          stddla[j].ion[Z].state[ion,idx].lambda = dumd1[k]
          stddla[j].ion[Z].state[ion,idx].sigclm = dumr2[k]
          stddla[j].ion[Z].state[ion,idx].flgclm = dumi1[k]
          stddla[j].ion[Z].state[ion,idx].flginst = dumi2[k]
      endfor

  endfor

return
end
