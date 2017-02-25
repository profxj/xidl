;+ 
; NAME:
; lls_fillion
;  V1.1
;
; PURPOSE:
;    Fills up the column densities for all of the ions observed for a
;    given LLS.  This program is unlikely to be called by anything
;    except lls_struct.
;
; CALLING SEQUENCE:
;   
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
;   18-May-2006 Written by JXP
;-
;------------------------------------------------------------------------------
pro lls_fillion, stdlls, nn, sysi, ROOT=root

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'lls_fillion, lls, nn, sysi, ROOT=(v1.1)' 
    return
  endif 

  if not keyword_set(ROOT) then root = ''
;

  close, 2
  dumc = ''
  dumi = 0

; Create the abundance array

  openr, 2, strtrim(root+stdlls[nn].systems[sysi].abndfil,2)
  readf, 2, dumc
  readf, 2, dumi
  stdlls[nn].systems[sysi].ndfil = dumi
  ;; HIRES
  if stdlls[nn].systems[sysi].ndfil MOD 2 GT 0 then begin
      readf, 2, dumc
      stdlls[nn].systems[sysi].Hfil = dumc
  endif
  ;; ESI 
  if stdlls[nn].systems[sysi].ndfil MOD 4 GT 1 then begin
      readf, 2, dumc
      stdlls[nn].systems[sysi].Efil = dumc
  endif
  ;; UVES
  if stdlls[nn].systems[sysi].ndfil MOD 8 GT 3 then begin
      readf, 2, dumc
      stdlls[nn].systems[sysi].Ufil = dumc
  endif
  ;; X?
  if stdlls[nn].systems[sysi].ndfil MOD 16 GT 7 then begin
      readf, 2, dumc
      stdlls[nn].systems[sysi].Xfil = dumc
  endif
  ;; MIKEB
  if stdlls[nn].systems[sysi].ndfil MOD 32 GT 15 then begin
      readf, 2, dumc
      stdlls[nn].systems[sysi].MBfil = dumc
  endif
  ;; MIKER
  if stdlls[nn].systems[sysi].ndfil MOD 64 GT 31 then begin
      readf, 2, dumc
      stdlls[nn].systems[sysi].MRfil = dumc
  endif
  readf, 2, dumc
  readf, 2, dumc
  close, 2

  ;; Allfil
  ionfil = root+strmid(dumc,0,strlen(strtrim(dumc,2))-3)+'all'
  if file_lines(ionfil) EQ 0 then return
  
  readcol, ionfil, format='I,I,F,F,I,I', dumi1, dumi2, dumr1, $
           dumr2, dumi3, dumi4, /silent

  ;; Central ion
  for k=0,n_elements(dumr1)-1 do begin
      stdlls[nn].systems[sysi].ion[dumi1[k]].state[dumi2[k]].clm = dumr1[k]
      stdlls[nn].systems[sysi].ion[dumi1[k]].state[dumi2[k]].sigclm = dumr2[k]
      stdlls[nn].systems[sysi].ion[dumi1[k]].state[dumi2[k]].flgclm = dumi3[k]
  endfor

  ;; Ionfil
  ionfil = root+strmid(dumc,0,strlen(strtrim(dumc,2))-3)+'ion'
  if file_lines(ionfil) EQ 0 then return

  readcol, ionfil, format='D,F,F,I,I', dumd1, dumr1, $
           dumr2, dumi1, dumi2, /silent
  
  ;; Central ion
  for k=0,n_elements(dumr1)-1 do begin
      ;; Grab atomic number and ion
      getion, dumd1[k], ion, Z=z
      stdlls[nn].systems[sysi].ion[Z].indx[ion]++
      idx = stdlls[nn].systems[sysi].ion[Z].indx[ion]
      stdlls[nn].systems[sysi].ion[Z].state[ion,idx].clm = dumr1[k]
      stdlls[nn].systems[sysi].ion[Z].state[ion,idx].lambda = dumd1[k]
      stdlls[nn].systems[sysi].ion[Z].state[ion,idx].sigclm = dumr2[k]
      stdlls[nn].systems[sysi].ion[Z].state[ion,idx].flgclm = dumi1[k]
      stdlls[nn].systems[sysi].ion[Z].state[ion,idx].flginst = dumi2[k]
  endfor

  return
end
