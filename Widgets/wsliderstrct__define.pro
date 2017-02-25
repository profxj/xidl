;+ 
; NAME:
; wsliderstrct__define
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
pro wsliderstrct__define

;  This routine defines the structure for the Slider widget

  device, get_screen_size=ssz
  fonts = x_widget_setfont(ssz[0])

  tmp = {wsliderstrct, $
         name: '',$
         font: '', $
         uname: '', $
         xsize: 0., $
         nslide: 0, $
         curslide: 0, $
         value: dblarr(30), $
         droplist: strarr(30), $
         min: dblarr(30), $
         max: dblarr(30), $
         fonts: fonts, $
         drop_id: 0L, $
         slider_id: 0L, $
         max_id:0L, $
         min_id: 0L, $
         base_id: 0L $
        }

end
  
         
