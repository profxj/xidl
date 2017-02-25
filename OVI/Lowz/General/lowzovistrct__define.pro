;+ 
; NAME:
; lowzovidat__define
;  V1.1
;
; PURPOSE:
;    Structure summarizing something (not sure what)
;
; CALLING SEQUENCE:
;   tmp = {lowzovistrct}
;
; INPUTS:
;
; RETURNS:
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
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Oct-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro lowzovistrct__define

  tmp = {lowzovistrct, $
         galstrct: ' ', $    ; Structure for galaxy info
         imgstrct: ' ', $    ; Structure for imaging info
         lyastrct: ' ', $    ; Lya absorbers
         ovistrct: ' ', $    ; OVI absorbers
         complete: 0., $     ; Completeness
         flg_complete: 0, $  ; Flag of completeness
         zqso: 0. $          ; QSO redshift
         }

end
  
