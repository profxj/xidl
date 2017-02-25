;+ 
; NAME:
; arclinstrct__define
;    Version 1.1
;
; PURPOSE:
;    Arc line structure
;
; CALLING SEQUENCE:
;   tmp = {arclinstrct}
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
;   Written by JXP
;-
;------------------------------------------------------------------------------
pro arclinstrct__define

;  This routine defines the exraction structure

  tmp = {arclinstrct, $
         wave: 0.d, $    ; Lab measurement
         flg_qual: 0, $    ; Quality: 0=NG, 1=Poor, 2=Fair, 3=G, 4=VG, 5= Exc
         flg_plt: 0, $     ; 1=Plot (and fit)
         pix: 0., $       ; Pixel value
         wave_fit: 0.d , $ ; Fit value
         name: ' ' $       ; Name
        }

end
  
         
