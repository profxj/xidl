;+ 
; NAME:
; x_arclist   
;    Version 1.0
;
; PURPOSE:
;    Reads a line list into an arclin strutcture
;
; CALLING SEQUENCE:
;   
;   x_arclist, linelist, lines
;
; INPUTS:
;   linelist  - Name of line list
;
; RETURNS:
;
; OUTPUTS:
;   arclinstr   -  Arc line structure
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_arclist, linelist, lines
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_arclist, linelist, lines, GDONLY=gdonly


;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_arclist, linelist, arclinstr, /GDONLY [v1.0]'
    return
  endif 

; Optional Keywords

  tmp = { arclinstrct }
  readcol, linelist, FORMAT='D,I,A', wav, flg, nam, /silent

  ; lines Structure
  nlin = n_elements(wav)
  lines = replicate(tmp, nlin)
  lines.wave = temporary(wav)
  lines.flg_qual = temporary(flg)
  lines.name = temporary(nam)

  if keyword_set( GDONLY ) then begin
      gdlin = where(lines.flg_qual NE 0)
      lines = lines[gdlin]
  endif

  return
end

