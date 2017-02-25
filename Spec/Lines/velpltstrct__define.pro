;+ 
; NAME:
; velpltstrct__define
;   Version 1.1
;
; PURPOSE:
;  Structure that is useful for velocity plots
;
; CALLING SEQUENCE:
;   tmp = {velpltstrct}
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
;---------------------------------------------------
pro velpltstrct__define

;  

  tmp = {velpltstrct, $
         wrest: 0.d, $
         ymnx: fltarr(2), $
         name: ' ', $
         flg: 0, $  ; 1=Plot
         blnd: fltarr(50,2) $
         }

end
  
         
