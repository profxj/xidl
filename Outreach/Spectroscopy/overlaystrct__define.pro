;+ 
; NAME:
; overlaystrct__define
;    Version 1.0
;
; PURPOSE:
;    Defines the structure for the Slider widget
;
; CALLING SEQUENCE:
;  tmp = {dlastruct}
;
; INPUTS:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Mar-2008 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro overlaystrct__define

;  This routine defines the structure for the Slider widget


  tmp = {overlaystrct, $
         name: '',$
         wave: dblarr(10000L), $
         npix: 0L, $
         fx: fltarr(10000L), $
         f_norm: 0., $
         param: fltarr(10) $  ;;  Ways to modify the overlay, e.g. T
        }

end
  
         
