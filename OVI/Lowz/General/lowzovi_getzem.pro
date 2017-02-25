;+ 
; NAME:
; lowzovi_getzem   
;   Version 1.0
;
; PURPOSE:
;    Returns zem from the structure using the partial name
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
function lowzovi_getzem, struct, QSO_NM

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
          'kcut = lowzovi_getzem(struct, QSO_NM) [v1.0]'
    return, -1
 endif 
  

  ;; INIT
  nm_len = strlen(QSO_NM)

  mt = where(strmatch(strmid(struct.qso,0,nm_len), QSO_NM), nmt)
  if nmt NE 1 then stop

  return, struct[mt].qso_zem
end

