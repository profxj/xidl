;+ 
; NAME:
; newabslinstrct__define
;   Version 1.1
;
; PURPOSE:
;  Current structure for absorption lines
;
; CALLING SEQUENCE:
;   tmp = {newabslinstrct}
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
pro nofabslinstrct__define

;  This routine is *only* for creating the linelst FITS file using
;  x_mkabslin.pro

  tmp = {nofabslinstrct, $
         ion: ' ', $
         wrest: 0.d, $
         fval: 0.d, $
         gamma: 0., $
         A: 0., $           ; Einstein coeff
         j: 0., $           ; Total ang mom (z projection of lower level)
         E: 0., $           ; Excitation energy (cm^-1)
         set: 0, $          ; Groups abs lines together
         N: 0., $           ; Colm (log) 
         Nsig: 0., $
         b: 0., $           ; Doppler parameter
         bsig: 0., $        ; Doppler parameter error
         zabs: 0.d,  $       ; Redshift
         zsig: 0.d   $       ; Error
         }

end
  
         
