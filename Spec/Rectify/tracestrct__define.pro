;+ 
; NAME:
; tracestrct__define
;   Version 1.1
;
; PURPOSE:
;  This routine defines the trace structure
;  Not in use
;
; CALLING SEQUENCE:
;   tmp = {tracestrct}
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
;----------------------------------------------------------------------
pro tracestrct__define

;  This routine defines the trace structure
;  Not in use

  tmp = {tracestrct, $
         xcen: ,  $       ; Name
         nord: 0,   $      ; Order number
         nrm: dblarr(2),   $      ; Normalization for x dimension
         lsig: 0., $
         hsig: 0., $
         niter: 0L, $
         minpt: 0L, $
         maxrej: 0L, $
         flg_rej: 0, $
         rms: 0., $
         ffit: ptr_new(/allocate_heap) $  ; 
         }

end
  
         
