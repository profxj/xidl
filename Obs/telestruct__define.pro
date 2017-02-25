; NAME:
; telestruct__define
;    Version 1.1
;
; PURPOSE:
;    Defines a structure for telescopes
;
; CALLING SEQUENCE:
;  tmp = {dlastruct}
;
; INPUTS:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Oct-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro telestruct__define

;  This routine defines the structure for direct images

  tmp = {telestruct, $
         name: '', $  ; KeckI, KeckII, Lick-3m 
         area: 0.,$
         plate_scale: 0. $
        }

end
  
         
