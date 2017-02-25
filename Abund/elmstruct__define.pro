;+ 
; NAME:
; elmstruct__define
;    Version 1.1
;
; PURPOSE:
;    Defines the structure for Elemental info
;
; CALLING SEQUENCE:
;  tmp = {elmstruct}
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
pro elmstruct__define

;  This routine defines the structure for direct images

  tmp = {elmstruct, $
         flgclm: 0,$  ; 0=Nothing; 1=Measurement; 2=Lower limit; 3=Upper
         flginst: 0,$ ; 1=HIRES, 2=ESI, 64=Other
         clm: 0.0,$   ; Logarithmic values (usually)
         sigclm: 0.0 }

end
  
         
