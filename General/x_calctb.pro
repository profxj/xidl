;+ 
; NAME:
; x_calctb   
;   Version 1.0
;
; PURPOSE:
;    Launches a cw_field and grabs input from the user
;
; CALLING SEQUENCE:
;   
;   num = x_guinum(flg)
;
; INPUTS:
;   flg = 0: float, 1: double, 2: Long
;
; RETURNS:
;   num - Number
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   num = x_guinum( 0 )
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_calctb, b, T, MA=MA, LOGT=logt

;common x_slctline_ans

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_calctb, b, t, MA='
    return
  endif 

;  Optional Keywords
  if not keyword_set( MA ) then MA = 1.

;; Constants
  c = x_constants(/cgs)

  T = (b*1e5)^2 * MA * c.mp / c.k
  if keyword_set( LOGT ) then T = alog10(T)

  return
;    
end
