;+ 
; NAME:
; printcol
;   Version 1.1
;
; PURPOSE:
;    Prints a series of arrays to the screen in column format.
;   The program will print as many entries as found in v1.
;
; CALLING SEQUENCE:
;   
;   printcol, v1, v2, [v3-v8], FORMAT=''
;
; INPUTS:
;   v1       - Vector 1
;   v2       - Vector 2
;   [v3-v8]  - Additional vectors.
;
; RETURNS:
;
; OUTPUTS:
;   Prints v1, v2 to screen
;
; OPTIONAL KEYWORDS:
;   FORMAT -  FORTRAN formatting
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The program keys off the number of elements in v1
;
; EXAMPLES:
;   printcol, array1, array2
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-June-2001 Written by JXP
;-
;------------------------------------------------------------------------------
pro printcol, v1, v2, v3, v4, v5, v6, v7, v8, FORMAT= format

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'printcol, v1, v2, [v3, v4, v5, v6], FORMAT= [v1.1]'
    return
  endif 

;

  flgvn = N_params()

  ;
  for i=0L,n_elements(v1)-1 do begin
      ;; Brute force
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
