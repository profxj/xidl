	;+ 
; NAME:
; printcol
;   Version 1.1
;
; PURPOSE:
;    Prints a series of arrays to the screen
;
; CALLING SEQUENCE:
;   
;   printcol, v1, v2, [v3-v6], FORMAT=''
;
; INPUTS:
;   v1       - Vector 1
;   v2       - Vector 2
;
; RETURNS:
;
; OUTPUTS:
;   Prints v1, v2 to screen
;
; OPTIONAL KEYWORDS:
;   v3-v6       - Vectors 3-6
;
; OPTIONAL OUTPUTS:
;   FORMAT -  FORTRAN formatting
;
; COMMENTS:
;   The program keys off the number of elements in v1
;
; EXAMPLES:
;   printcol, array1, array2
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-June-2001 Written by JXP
;-
;------------------------------------------------------------------------------
pro printcol, v1, v2, v3, v4, v5, v6, v7, v8, FORMAT= format


; writecol -- Writes a 2 column ascii file

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'printcol, v1, v2, [v3, v4, v5, v6] FORMAT= '
    return
  endif 

;

  flgvn = 2
  if keyword_set( v8 ) then    flgvn    = flgvn + 1
  if keyword_set( v7 ) then    flgvn    = flgvn + 1
  if keyword_set( v6 ) then    flgvn    = flgvn + 1
  if keyword_set( v5 ) then    flgvn    = flgvn + 1
  if keyword_set( v4 ) then    flgvn    = flgvn + 1
  if keyword_set( v3 ) then    flgvn    = flgvn + 1

;

  for i=0,n_elements(v1)-1 do begin
      case flgvn of 
          8: print, FORMAT = format, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i]
          7: print, FORMAT = format, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i]
          6: print, FORMAT = format, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i]
          5: print, FORMAT = format, v1[i], v2[i], v3[i], v4[i], v5[i]
          4: print, FORMAT = format, v1[i], v2[i], v3[i], v4[i] 
          3: print, FORMAT = format, v1[i], v2[i], v3[i]
          2: print, FORMAT = format, v1[i], v2[i]
          else: stop
      endcase
  endfor

return
end
