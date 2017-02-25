;+ 
; NAME:
; ionstruct__define
;    Version 1.1
;
; PURPOSE:
;    Defines the structure for ion abundance info
;
; CALLING SEQUENCE:
;  tmp = {ionstruct}
;
; INPUTS:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   1-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ionstruct__define


  tmp1 = { ionsubstr, $
           flgclm: 0,$          ; 0=Nothing; 1=Measurement; 2=Lower; 3=Upper
           flginst: 0,$         ; 1=HIRES, 2=ESI, 64=Other
           lambda: 0.d, $
           vmnx: fltarr(2), $ ; Limits for integration
           clm: 0.0,$           ; Logarithmic
           sigclm: 0.0 }

  tmp = {ionstruct, $
         indx: lonarr(16), $
         state: replicate(tmp1, 61, 30) $  ; I*** in 31.., II*** in 41..
         }

end
  
         
