;+ 
; NAME:
; x_qabnd
;   Version 1.0
;
; PURPOSE:
;    Finds data point of an array nearest the input value
;      and returns the array member
;
; CALLING SEQUENCE:
;   
;   indx = x_qabnd(val, xdat, XVAL=)
;
; INPUTS:
;   xdat - Data
;   val  - value
;
; RETURNS:
;   indx  - Index in the array
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   XVAL - Value of the array at that index
;
; COMMENTS:
;
; EXAMPLES:
;   indx = x_qabnd( 1.0, findgen(1000))
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   19-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_qabnd, file, NHI

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'x_qabnd, file, [NHI]  [V1.0]'
    return
  endif 

  if not keyword_set(NHI) then NHI = 21.0
;; Parse the file

  readcol, file, Z, NX, sigN, flg, dum, FORMAT='I,F,F,I,I'

  ;; LOOP
  for i=0L,n_elements(Z)-1 do begin
      ;; Get abnd
      getabnd, elm, Z[i], abnd, /flag

      ;; Calculate [X/H]
      XH = NX[i] - NHI + 12. - abnd

      ;; PRINT
      case flg[i] of 
          1: print, elm, XH, FORMAT='(a2,3x,f6.3)'
          2: print, elm, '>',XH, FORMAT='(a2,1x,a1,1x,f6.3)'
          3: print, elm, '<',XH, FORMAT='(a2,1x,a1,1x,f6.3)'
          else: stop
      endcase
  endfor
  return
end

