;+ 
; NAME:
; stdstruct__define
;   Version 1.1
;
; PURPOSE:
;   Structure to handle observations of Standard stars
;
; CALLING SEQUENCE:
;  tmp = {stdstruct}
;
; INPUTS:
;
; RETURNS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro stdstruct__define

;  This routine defines the structure for standard stars

  tmp = {stdstruct, $
         flg_anly: 0,  $    ; Flag to include in analysis (0=No)
         Name: ' ',     $    ; Maps image onto entire window  
         filter: ' ',   $    ; Filter (string)
         Mag: 0.,       $    ; Magnitude
         sig_Mag: 0.,   $    ; Sigma in Magnitude
         AM: 0.,        $    ; Air Mass
         CLR: 0.,       $    ; Color Term (Value)
         sCLR: ' '      $
         }

end
  
         
