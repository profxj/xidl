;+ 
; NAME:
; lowzovi_kstarcut   
;   Version 1.0
;
; PURPOSE:
;    Launches a cw_field and grabs input from the user
;
; CALLING SEQUENCE:
;   
;   string = lowzovi_kstarcut(title)
;
; INPUTS:
;   title - Title
;
; RETURNS:
;   string - String
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
;   area = lowzovi_kstarcut( rmag, /INIT)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function lowzovi_kstarcut, rmag, INIT=init

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'kcut = lowzovi_kstarcut(rmag, /INIT) [v1.0]'
    return, -1
  endif 

; INIT
