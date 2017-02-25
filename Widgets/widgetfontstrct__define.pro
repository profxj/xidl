;+ 
; NAME:
; widgetfontstrct__define
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
pro widgetfontstrct__define

;  This routine defines the structure for the Slider widget


  tmp = {widgetfontstrct, $
         scr_xpix: 0L, $
         tiny_font: '', $
         small_font: '', $
         big_font: '', $
         huge_font: '' $
        }

end
  
         
