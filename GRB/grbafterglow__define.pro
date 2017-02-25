;+ 
; NAME:
; grbafterglow__define
;   Version 1.1
;
; PURPOSE:
;  Structure for a GRB light curve
;
; CALLING SEQUENCE:
;   tmp = {abslinstrct}
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
pro grbafterglow__define

;  This routine defines the line list structure

  tmp = {grbafterglow, $
         nam: ' ', $
         ref: strarr(10), $
         refnum: lonarr(10), $
         z: 0., $
         alpha: 0., $
         beta: 0., $
         t0: 0., $
         mag: 0., $             ; Magnitude at t0
         mag_wv: 0.  $          ; Wavelength of the filter
         }

end
  
         
