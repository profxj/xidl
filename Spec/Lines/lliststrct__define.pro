;+ 
; NAME:
; lliststrct__define
;   Version 1.1
;
; PURPOSE:
;  Structure for the line lists
;
; CALLING SEQUENCE:
;   tmp = {lliststrct}
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
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;  Written by JXP
;-
;----------------------------------------------------------------------
pro lliststrct__define

;  This routine defines the line list structure

  tmp = {lliststrct, $
         name: ' ',  $       ; Name
         elm: ' ',  $        ; Elm
         ion: ' ',  $        ; Ion
         wave: 0.d,   $      ; Wavelength
         fval: 0.d,       $  ; f-value   
         gamma: 0.,       $  ; gamma
         ref: 0, $           ; Reference
         flg: 0  $           ; Flag
         }

end
  
         
