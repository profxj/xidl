;+ 
; NAME:
; fill_elmxh   
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
pro fill_elmxh, stddla, NOHIS=nohis, ROOT=root

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'fill_elmxh, dla, /NOHIS, ROOT= [v1.1]'
    return
  endif 

  if not keyword_set(ROOT) then root = ''

  close, 2
  dumc = ''
  dumi = 0

; Create the abundance array

  ndla = n_elements(stddla)

  for j=0,ndla-1 do begin

      ;; HI
      stddla[j].XH[1].flgclm = 1
      stddla[j].XH[1].flginst = 1
      stddla[j].XH[1].clm = 0.
      if stddla[j].sigNHI[1] GT 0. then $
        stddla[j].XH[1].sigclm = mean(stddla[j].sigNHI) else $
        stddla[j].XH[1].sigclm = stddla[j].sigNHI[0]

      ;; Abund fil
      a = findfile(root+strtrim(stddla[j].abndfil,2),count=na)
      if na EQ 0 then begin
          print, 'fill_elmxh: No Abund file -- ', root+stddla[j].abndfil
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
      XHfil = strmid(dumc,0,strlen(dumc)-3)+'XH'
      
      a = file_search(root+XHfil,count=na)
      if na EQ 0 then begin
          close, 2
          continue
      endif
      readfmt, root+XHfil, 'i2,1x,f6.3,2x,f5.3,1x,i1,2x,i2',$
        dumi1, dumr1, dumr2, dumi2, dumi3, /silent
      for k=0,n_elements(dumi1)-1 do begin
          ;; Elm
          stddla[j].elm[dumi1[k]].clm = dumr1[k] 
          stddla[j].elm[dumi1[k]].sigclm = dumr2[k]
          stddla[j].elm[dumi1[k]].flgclm = dumi2[k]
          ;; Instrument
          stddla[j].elm[dumi1[k]].flginst = dumi3[k]
          stddla[j].XH[dumi1[k]].flginst = dumi3[k]
          ;; XH
          getabnd, nm, dumi1[k], abnd, flag=1
          stddla[j].XH[dumi1[k]].clm = dumr1[k] - stddla[j].NHI + 12.0 - abnd
          ;; Sigma
          if not keyword_set( NOHIS ) then begin 
              if stddla[j].sigNHI[1] GT 0. then $
                sigNHI = mean(stddla[j].sigNHI) else $
                sigNHI = stddla[j].sigNHI[0]
              stddla[j].XH[dumi1[k]].sigclm = sqrt(dumr2[k]^2 + sigNHI^2)
          endif else stddla[j].XH[dumi1[k]].sigclm = dumr2[k]
          ;; Flag
          stddla[j].XH[dumi1[k]].flgclm = dumi2[k]
      endfor
      close, 2
  endfor

return
end
