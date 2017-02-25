;+ 
; NAME:
; abslinstrct__define
;   Version 1.1
;
; PURPOSE:
;  Structure for a simple absorption line. 
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
pro abslinstrct__define

;  This routine defines the line list structure

  tmp = {abslinstrct, $
         ion: ' ', $
         wrest: 0.d, $
         f: 0.d, $
         gamma: 0., $
         set: 0, $          ; Groups abs lines together
         N: 0., $           ; Colm (log) 
         Nsig: 0., $
         b: 0., $
         bsig: 0., $
         zabs: 0.d,  $       ; Redshift
         zsig: 0.d  $       ; Error
         }

end
  
         
