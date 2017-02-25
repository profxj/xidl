;+ 
; NAME:
; getion
;  Version 1.1
;
; PURPOSE:
;  Given an atomic wavelength, return the ionic value
;    (e.g. 1, 2, 3).  Uses getfnam to grab the transition name.  
;
; CALLING SEQUENCE:
;   
;   getion, wave, ion, [elm], Z=z, ABND=abnd, NM=nm, INM=inm, FNM=fnm
;
; INPUTS:
;   wave       - rest wavelength of transition
;
; RETURNS:
;
; OUTPUTS:
;   ion        - Ion value (e.g. 1,2,3)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   [elm]   - Name of element 
;   Z=      - Atomic number of the element
;   ABND=   - Atomic abundance on the 12 log scale
;   NM=     - Name of the ion (e.g. IV)
;   FNM=    - Full name of the ion (e.g. CIV)
;   /INM    - Indicates wave is actually the Name of the transition
;
; COMMENTS:
;
; EXAMPLES:
;   getion, 1215.6701, ion, elm
;
; PROCEDURES CALLED:
;  getfnam
;
; REVISION HISTORY:
;   11-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro getion, wave, ion, elm, Z=z, ABND=abnd, NM=nm, INM=inm, FNM=fnm

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'getion, wave, ion, [elm], Z=, ABND=, NM=, /INM, FNM= [v1.1]'
    return
  endif 

;
  if not keyword_set(INM) then getfnam, wave, fval, nam $
  else nam = wave

;     Parse nam
  for i=0,10 do begin
      if strmid(nam,i,1) NE ' '  then begin
          c1 = i
          break
      endif   
  endfor

  if strmid(nam,c1+1,1) NE 'I' AND strmid(nam,c1+1,1) NE 'V' then c2 = c1+2 $
  else c2 = c1+1

  elm = strmid(nam,c1,c2-c1)

; ION code
  cion = strmid(nam,c2)
  ipos = strpos(cion, ' ')
  if ipos EQ -1 then ipos = strlen(cion)
  nm = strmid(cion, 0, ipos)
;  nm = strtrim(strmid(nam,c1+1,3),2)
  case strtrim(nm,2) of 
      'I': ion = 1
      'I*': ion = 1
      'I**': ion = 1
      'II': ion = 2
      'II*': ion = 2
      'II**': ion = 2
      'III': ion = 3
      'IV': ion = 4
      'V': ion = 5
      'VI': ion = 6
      'VII': ion = 7
      'VIII': ion = 8
      'IX': ion = 9
      'X': ion = 10
      'XI': ion = 11
      'XII': ion = 12
      else: stop
  endcase

  ;; Z value
  if arg_present(Z) or arg_present(ABND) or arg_present(FNM) then begin
      getabnd, elm, Z, abnd
      fnm = elm+nm
  endif

  return
end
