;+ 
; NAME:
; x_convarclst   
;    Version 1.1
;
; PURPOSE:
;    Converts an arc line list (e.g. Murphy) to my format
;
; CALLING SEQUENCE:
;   x_convarclst, infil, outfil, FLG=
;
; INPUTS:
;   infil -- Input line list
;
; RETURNS:
;
; OUTPUTS:
;   outfil -- Name of output file
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  arclinstrct__define
;  readcol
;
; REVISION HISTORY:
;   17-Sep-2007 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_convarclst, infil, outfil, FLG=flg


;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_convarclst, infil, outfil, FLG= [v1.1]'
    return
  endif 

; Optional Keywords
  if not keyword_set(FLG) then flg = 1

  ;; Read
  case flg of 
      1: begin
          readcol, infil, dumf, wav, dumf, nam, format='D,D,F,A'
          flg = replicate(1, n_elements(wav))
      end
      else: stop
  endcase

  writecol, outfil, wav, flg, nam, FMT='(f11.6,1x,i1,1x,a4)'

  return
end

