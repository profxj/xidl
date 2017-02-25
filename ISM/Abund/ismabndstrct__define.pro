;+ 
; NAME:
; ismabndstrct__define
;   Version 1.1
;
; PURPOSE:
;  Structure for COG analysis
;
; CALLING SEQUENCE:
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
pro ismabndstrct__define

;  This routine defines a set of parameters for calculating the
;  column/EW of a given absorption line

  tmp = {ismabndstrct, $
         datfil: ' ', $
         ionnm: ' ', $
         Z: 0, $
         ion: 0, $
         wrest: 0.d, $
         zabs: 0.d, $
         dv: dblarr(2), $  ; Velocity interval for the calculation
         fval: 0.d, $    ; f value
         flg_anly: 0, $ ; Analysis flag (0=No, 1=Measure but dont use, 2=Full)
         flg_eye: 0, $ ; Eyeball flag (0=No comment, 1=detected, 2=blended, 4=bad continuum) [Binary] 
         flg_limit: 0, $ ; Limit flag (1=Normal, 2=Lower, 3=Upper)
         inst: 0, $ ; Instrument flag
         AODM: 0.d, $ ;  AODM column
         sigAODM: 0.d, $ ; 
         boxEW: 0., $  ; Box car EW
         sigboxEW: 0., $  ; Error
         gauss_fit: dblarr(3), $ ; Gaussian fit
         gauss_sig: dblarr(3), $ ; Gaussian error
         EW: 0.d, $ ; Ang (rest)
         sigEW: 0.d  $ ; Ang (rest)
         }

end
  
         
