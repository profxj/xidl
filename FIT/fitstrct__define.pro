;+ 
; NAME:
; fitstrct__define
;   Version 1.1
;
; PURPOSE:
;    IDL structure for 1D fitting
;
; CALLING SEQUENCE:
;   
;   tmp = {fitstrct}
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
;   tmp = {fitstrct}
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   31-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fitstrct__define

;  This routine defines the Fit structure
  tmp = {fitstrct, $
         func: ' ',  $       ; Name
         nord: 0L,   $      ; Order number
         nrm: dblarr(2),   $      ; Normalization for x dimension
         lsig: 0., $
         hsig: 0., $
         niter: 0L, $
         minpt: 0L, $
         maxrej: 0L, $
         flg_rej: 0, $
         rms: 0.d, $
         ffit: ptr_new(/allocate_heap) $  ; 
         }

end
  
         
