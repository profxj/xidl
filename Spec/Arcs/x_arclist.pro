;+ 
; NAME:
; x_arclist   
;    Version 1.1
;
; PURPOSE:
;    Reads a line list into an arclin strutcture
;
; CALLING SEQUENCE:
;   x_arclist, linelist, lines, /GDONLY
;
; INPUTS:
;   linelist  - Name of line list (formatted: D, I, A)
;
; RETURNS:
;   arclinstr   -  Arc line structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /GDONLY -- Restrict the line list to those where flg_qual NE 0
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_arclist, linelist, lines
;
; PROCEDURES/FUNCTIONS CALLED:
;  arclinstrct__define
;  readcol
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
             'x_arclist, linelist, arclinstr, /GDONLY [v1.1]'
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

