;+ 
; NAME:
; cogstrct__define
;   Version 1.1
;
; PURPOSE:
;  Structure for COG analysis
;
; CALLING SEQUENCE:
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
;------------------------------------------------------------------------------
pro cogstrct__define

;  This routine defines the line list structure

  tmp = {cogstrct, $
         ionnm: ' ', $
         Z: 0, $
         ion: 0, $
         nlin: 0L, $
         wrest: replicate(1.d,30), $  ; Ang (rest)
         dv: dblarr(30,2), $  ; Velocity interval for the calculation
         f: dblarr(30), $    ; f value
         inst: intarr(30), $ ; Instrument flag
         EW: dblarr(30), $ ; Ang (rest)
         sigEW: dblarr(30) $ ; Ang (rest)
         }

end
  
         
