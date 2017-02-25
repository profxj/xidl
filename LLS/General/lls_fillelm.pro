;+ 
; NAME:
; lls_fillelm   
;   Version 1.1
;
; PURPOSE:
;  Given a DLA structure, fills up the X/H info.  Simply reads in the
;  info from the .XH files.  This program is not likely to be called
;  by any program except parse_dlalst.
;
; CALLING SEQUENCE:
;   
;  fill_elmxh, dla
;
; INPUTS:
;  dla  -- IDL DLA structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOHIS -- Do not include HI error in analysis
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   29-May-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lls_fillelm, stdlls, nn, sys, NOHIS=nohis, ROOT=root

  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'lls_fillelm, lls, nn, sys, /NOHIS, ROOT= [v1.1]'
    return
  endif 

  if not keyword_set(ROOT) then root = ''

  close, 2
  dumc = ''
  dumi = 0

; Create the abundance array

  ;; Abund fil
  a = findfile(strtrim(root+stdlls[nn].systems[sys].abndfil,2),count=na)
  if na EQ 0 then begin
      print, 'lls_fillelm: No Abund file -- ', root+stdlls[nn].systems[sys].abndfil
      return
  endif
  
  openr, 2, strtrim(root+stdlls[nn].systems[sys].abndfil,2)
  readf, 2, dumc
  readf, 2, dumi
  stdlls[nn].systems[sys].ndfil = dumi
  
  ;; HIRES
  if stdlls[nn].systems[sys].ndfil MOD 2 GT 0 then begin
      readf, 2, dumc
      stdlls[nn].systems[sys].Hfil = dumc
  endif
  ;; ESI 
  if stdlls[nn].systems[sys].ndfil MOD 4 GT 1 then begin
      readf, 2, dumc
      stdlls[nn].systems[sys].Efil = dumc
  endif
  ;; UVES
  if stdlls[nn].systems[sys].ndfil MOD 8 GT 3 then begin
      readf, 2, dumc
      stdlls[nn].systems[sys].Ufil = dumc
  endif
  ;; X?
  if stdlls[nn].systems[sys].ndfil MOD 16 GT 7 then begin
      readf, 2, dumc
      stdlls[nn].systems[sys].Xfil = dumc
  endif
  ;; MIKEB
  if stdlls[nn].systems[sys].ndfil MOD 32 GT 15 then begin
      readf, 2, dumc
      stdlls[nn].systems[sys].MBfil = dumc
  endif
  ;; MIKER
  if stdlls[nn].systems[sys].ndfil MOD 64 GT 31 then begin
      readf, 2, dumc
      stdlls[nn].systems[sys].MRfil = dumc
  endif
  readf, 2, dumc
  readf, 2, dumc
  close, 2
  
  
  xhfil = strmid(dumc,0,strlen(strtrim(dumc,2))-3)+'XH'
  if file_lines(root+XHfil) EQ 0 then return
  readfmt, root+XHfil, 'i2,1x,f6.3,1x,f6.3,1x,f6.3,1x,i1,1x,i3', $
           dumi1, dumr1, dumr2, dumr3, dumi2, dumi3, /silent

  ;; XH
  for k=0,n_elements(dumi1)-1 do begin
      stdlls[nn].systems[sys].XH[dumi1[k]].clm = dumr1[k] 
      stdlls[nn].systems[sys].XH[dumi1[k]].sigclm = [dumr2[k], dumr3[k]]
      stdlls[nn].systems[sys].XH[dumi1[k]].flgclm = dumi2[k]
      stdlls[nn].systems[sys].XH[dumi1[k]].flginst = dumi3[k]
  endfor

  ;; Fill elm with .all

return
end
