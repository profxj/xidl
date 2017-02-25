;+ 
; NAME:
; lls_fillew
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
pro lls_fillew, stdlls, nn, sysi, ROOT=root

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'lls_fillew, lls, nn, sysi [v1.0]' 
    return
  endif 

  if not keyword_set(ROOT) then root = ''
;

  close, 2
  dumc = ''
  dumi = 0

; Create the abundance array


  ;; EW fil
  a = findfile(strtrim(ROOT+stdlls[nn].systems[sysi].abndfil,2),count=na)
  if na EQ 0 then begin
      print, 'lls_fillew: No Abund file -- ', $
             ROOT+stdlls[nn].systems[sysi].abndfil
      return
  endif

  ;; Open
  openr, 2, ROOT+strtrim(stdlls[nn].systems[sysi].abndfil,2)
  readf, 2, dumc
  readf, 2, dumi
;  stdlls[nn].systems[sysi].ndfil = dumi
  if(dumi mod 2 GT 0) then begin
      readf, 2, dumc 
      stdlls[nn].systems[sysi].Hfil = dumc
  endif
  if(dumi mod 4 GT 1) then begin
      readf, 2, dumc 
      stdlls[nn].systems[sysi].Efil = dumc
  endif
  if(dumi mod 8 GT 3) then begin
      readf, 2, dumc 
      stdlls[nn].systems[sysi].Ufil = dumc
  endif
  if(dumi mod 16 GT 7) then begin
      readf, 2, dumc 
      stdlls[nn].systems[sysi].Xfil = dumc
  endif
  if(dumi mod 64 GT 31) then begin
      readf, 2, dumc 
      stdlls[nn].systems[sysi].Xfil = dumc
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
  
  ;; EW fil
  ionfil = root+strmid(dumc,0,strlen(strtrim(dumc,2))-3)+'EW'
  a = findfile(ionfil,count=na)
  if na EQ 0 then begin
      print, 'lls_fillew: No EW file -- ', ionfil
      return
  endif
  if file_lines(ionfil) EQ 0 then return
  
  readcol, ionfil, format='D,F,F,I,I', wrest, EW, sigEW, inst, /sile
  
  ;; Central ion
  for k=0,n_elements(wrest)-1 do begin
      ;; Grab atomic number and ion
      getion, wrest[k], ion, Z=z
      stdlls[nn].systems[sysi].ion[Z].indx[ion]++
      idx = stdlls[nn].systems[sysi].ion[Z].indx[ion]
      stdlls[nn].systems[sysi].ion[Z].state[ion,idx].clm = EW[k] ; Rest EW
      stdlls[nn].systems[sysi].ion[Z].state[ion,idx].lambda = wrest[k]
      stdlls[nn].systems[sysi].ion[Z].state[ion,idx].sigclm = sigEW[k]
      stdlls[nn].systems[sysi].ion[Z].state[ion,idx].flginst = inst[k]
  endfor

return
end
