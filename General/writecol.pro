;+ 
; NAME:
; writecol
;
; PURPOSE:
;    Prints a series of arrays to a file in ASCII format
;
; CALLING SEQUENCE:
;   
;   writecol, filename, v1, v2, FMT=''
;
; INPUTS:
;   file     - Name of the ASCII file
;   v1       - Vector 1
;   v2       - Vector 2
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
;   v3-v9       - Vectors 3-9
;
; COMMENTS:
;   The program keys off the number of elements in v1
;
; EXAMPLES:
;   writecol, 'arrays.dat', array1, array2
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-June-2001 Written by JXP
;-
;------------------------------------------------------------------------------
pro writecol, file, v1, v2, v3, v4, v5, v6, v7, v8, v9, $
              v10, v11, v12, v13, v14, FMT=fmt, FILNUM=filnum


; writecol -- Writes a 2 column ascii file

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'writecol, file, v1, v2, [v3, v4, v5, v6, v7, v8, v9, v10, v11, ] FMT= '
    return
  endif 

;

  flgvn = 2
  if keyword_set( v14 ) then    flgvn    = flgvn + 1
  if keyword_set( v13 ) then    flgvn    = flgvn + 1
  if keyword_set( v12 ) then    flgvn    = flgvn + 1
  if keyword_set( v11 ) then    flgvn    = flgvn + 1
  if keyword_set( v10 ) then    flgvn    = flgvn + 1
  if keyword_set( v9 ) then    flgvn    = flgvn + 1
  if keyword_set( v8 ) then    flgvn    = flgvn + 1
  if keyword_set( v7 ) then    flgvn    = flgvn + 1
  if keyword_set( v6 ) then    flgvn    = flgvn + 1
  if keyword_set( v5 ) then    flgvn    = flgvn + 1
  if keyword_set( v4 ) then    flgvn    = flgvn + 1
  if keyword_set( v3 ) then    flgvn    = flgvn + 1
  if not keyword_set( FMT ) then    flgfmt    = 0 else begin
      flgfmt = 1 
      fmt = fmt[0]
  endelse

;

  if not keyword_set(FILNUM) then begin
      filnum = 1
      close, filnum
      openw, filnum, file
      flg_fil = 1
  endif

  for i=0L,n_elements(v1)-1 do begin
      case flgvn of 
          14: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i],v14[i]
          13: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i],v13[i]
          12: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i],v12[i]
          11: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i],v11[i]
          10: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i],v10[i]
          9: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],$
            v8[i],v9[i]
          8: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i]
          7: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i]
          6: printf, filnum, FORMAT=fmt, v1[i],v2[i],v3[i],v4[i],v5[i],v6[i]
          5: printf, filnum, FORMAT=fmt, v1[i], v2[i], v3[i], v4[i], v5[i]
          4: printf, filnum, FORMAT=fmt, v1[i], v2[i], v3[i], v4[i] 
          3: printf, filnum, FORMAT=fmt, v1[i], v2[i], v3[i]
          2: printf, filnum, FORMAT=fmt, v1[i], v2[i]
          else: stop
      endcase
  endfor
  if keyword_set(FLG_FIL) then close, filnum


return
end
