
; NAME:
; instrstruct__define
;    Version 1.1
;
; PURPOSE:
;    Defines a structure for an instrument (spectrometer)
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
pro instrstruct__define

;  This routine defines the structure for direct images

  tmp = {instrstruct, $
         name: '', $
         mag_perp: 0.,$
         mag_para: 0.,$
         pixel_size: 0.,$
         scale_perp: 0.,$
         scale_para: 0.0,$
         R: 0.,$
         mlambda: 0.0,$
         dark: 0.0,$
         readno: 0.0,$
         dichroic: '',$
         grating: '',$
         cwave: 0., $    ; Central wavelength
         swidth: 0.,$    ; Slit width
         sheight: 0.,$   ; Slit height
         wvmnx: fltarr(2),$   ; Slit height
         bins: 0,$       ; Spatial binning
         bind: 0,$       ; Dispersion binning
         dely: 0.0 $
        }

end
  
         
