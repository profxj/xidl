;+ 
; NAME:
; getion
;
; PURPOSE:
;    Given an atomic wavelength, return the fvalue and name
;
; CALLING SEQUENCE:
;   
;   getfnam, wave, fval, [nam]
;
; INPUTS:
;   wave       - ionic transition
;
; RETURNS:
;   f          - oscillator strength
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   nam        - Name of transition
;
; COMMENTS:
;
; EXAMPLES:
;   getfnam, 1215.6701, fval, name
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   11-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro getion, wave, ion, elm, Z=z

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'getfnam, wave, ion, [elm]'
    return
  endif 

;
  getfnam, wave, fval, nam

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
  nm = strmid(cion, 0, ipos)
;  nm = strtrim(strmid(nam,c1+1,3),2)
  case nm of 
      'I': ion = 1
      'II': ion = 2
      'III': ion = 3
      'IV': ion = 4
      'V': ion = 5
      'VI': ion = 6
      'VII': ion = 7
      'VIII': ion = 8
      else: stop
  endcase

  ;; Z value
  if arg_present(Z) then begin
      getabnd, nam, Z, abnd
  endif

  return
end
