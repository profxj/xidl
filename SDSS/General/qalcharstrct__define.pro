;+ 
; NAME:
; qalcharstrct__define
;   Version 1.1
;
; PURPOSE:
;  Defines the QAL structure for SDSS
; 
; CALLING SEQUENCE:
; tmp = {qalcharstrct}
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
;   30-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro qalcharstrct__define

tmp={qalcharstrct, $
     file_name: ' ', $
     qso_name: ' ', $
     plate: 0L, $
     fiberid: 0L, $
     qso_mag: 0.d, $ 
     qso_flux: 0.d, $   ;; 1300-1350 in the QSO rest frame (units 1e17 flamb)
     RA: 0.d, $
     DEC: 0.d, $
     Z_qso: 0., $
     SNR: 0.d, $
     flg_bal: 0, $   ; 1=BAL; 2=BAL; 3=Not a QSO
     flg_LLS: 0, $
     LLS_zabs: 0.d, $
     start_wave: 0.d, $
     nDLA1: 0, $         ;  Number of HI detections
     DLA_zabs1: dblarr(100), $
     DLA_score1: fltarr(100), $
     nDLA2: 0, $         ;  Number of metal-line systems (not necess DLA)
     DLA_zabs2: dblarr(100), $
     DLA_score2: fltarr(100), $
     DLA_hits:  fltarr(100), $  ; Number of DLA candidates altogether
     DLA_z: dblarr(100), $
     DLA_quality: fltarr(100) $ ; Combined quality from HI and metals
     }

return
end
     
