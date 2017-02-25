;+ 
; NAME:
; lndltstr__define
;   Version 1.1
;
; PURPOSE:
;   Structure to handle Landolt star data
;
; CALLING SEQUENCE:
;  tmp = {lndltstr}
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
;   07-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lndltstr__define

;  This routine defines the structure for Landolt data

  tmp = {lndltstr, $
         Name: ' ',     $    ; Maps image onto entire window  
         RA: ' ',       $
         DEC: ' ',      $
         V: 0.,         $
         BV: 0.,        $
         UB: 0.,        $
         VR: 0.,        $
         RI: 0.,        $
         VI: 0.,        $
         nobs: 0,       $    ; Number of times observed
         mobs: 0,       $    ; Number of nights observed
         sig_V: 0.,     $
         sig_BV: 0.,    $
         sig_UB: 0.,    $
         sig_VR: 0.,    $
         sig_RI: 0.,    $
         sig_VI: 0.     $
         }

end
  
         
