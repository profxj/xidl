;+ 
; NAME:
; h2linstrct__define
;  (V1.1)
;
; PURPOSE:
;    Molecular hydrogen structure
;
; CALLING SEQUENCE:
;   
;   tmp = {h2linstrct}
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
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   09-Feb-2004 Written by JXP
;-
;------------------------------------------------------------------------------
pro h2linstrct__define

;  This routine defines the line list structure

  tmp = {h2linstrct, $
         wrest: 0.d, $
         f: 0.d, $
         gamma: 0., $
         el: 0, $
         np: 0, $
         npp: 0, $
         Jp: 0, $
         Jpp: 0, $
         label: '', $
         set: 0, $          ; Groups abs lines together
         N: 0., $           ; Colm (log) 
         Nsig: 0., $
         b: 0., $
         bsig: 0., $
         zabs: 0.d,  $       ; Redshift
         zsig: 0.d  $       ; Error
         }

end
  
         
