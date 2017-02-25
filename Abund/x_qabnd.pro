;+ 
; NAME:
; x_qabnd
;   Version 1.1
;
; PURPOSE:
;    Given an abundance file, quickly converts to standard [X/H] notation
;
; CALLING SEQUENCE:
;   
;   x_qabnd, fil, [NHI]
;
; INPUTS:
;   fil - Abundance file (extension .XH generally)
;   NHI -- NHI value of system
;
; RETURNS:
;
; OUTPUTS:
;  Prints the [X/H] values to the screen
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_qabnd, 'q0405.XH', NHI=15.4
;
; PROCEDURES CALLED:
;  getabnd
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

