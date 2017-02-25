;+ 
; NAME:
; fit2dstrct__define
;   Version 1.1
;
; PURPOSE:
;    IDL structure for 2D surface fitting
;
; CALLING SEQUENCE:
;   
;   tmp = {fit2dstrct}
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
;   tmp = {fit2dstrct}
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   31-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fit2dstrct__define

;  This routine defines the Fit structure for a 2D surface fit

  tmp = {fit2dstrct, $
         func: ' ',  $       ; Name
         nx: 0,   $      ; Order number
         ny: 0,   $      ; Order number
         nrm: dblarr(2,2),   $   ; Normalization for x,y dimensions
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
  
         
