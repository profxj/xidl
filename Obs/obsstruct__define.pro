
; NAME:
; obsstruct__define
;    Version 1.1
;
; PURPOSE:
;    Defines a structure for observing
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
pro obsstruct__define

;  

  tmp = {obsstruct, $
         seeing: 0.,$
         airmass: 0.,$
         mphase: 0,$
         mstar: 0.,$       ; Mag
         filter: '',$      ; Mag flag:  0=AB
         mtype: 0,$        ; Mag flag:  1=AB, 2=Johnson
         exptime: 0.0, $
         redshift: 0.0, $
         template: '', $
         vega_template: '' $
        }

end
  
         
