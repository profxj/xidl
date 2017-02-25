;+ 
; NAME:
; civcandstrct__define
;  V1.1
;
; PURPOSE:
;    Defines the structure for HST CIV candidate search
;
; CALLING SEQUENCE:
;  tmp = {civcandstrct}
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
;   Feb-2008   created by JXP, adapted from fuselinstrct__define
;   Oct-2010   Increased AODM arrays to 801 pixels
;-
;------------------------------------------------------------------------------
pro civcandstrct__define

;  This routine defines the line list structure

  tmp = {civcandstrct, $
         QSO: ' ', $  ;   Name of Quasar
         RA:  0.d, $  ;   RA of QSO  (decimal degrees)
         DEC: 0.d, $  ;   DEC of QSO (decimal degrees)
         zqso: 0.d,$  ;   Redshift of quasar
         search_fil: ' ',$     ; File with all detected abs lines
         instr_fil: ' ', $     ; List of spectra filenames
         flg_sys: intarr(2), $         ; Flag describing CIV 'quality'
         rating_eye: 0, $              ; Eyeball rating
         comment_eye: bytarr(20L), $   ; Button comments (see the GUI)
         other_comments: '', $         ; Written comment
         aodm_vel: fltarr(801,2), $    ; AODM velocity array
         aodm_civ: fltarr(801,2), $  ; AODM array for CIV
         sigaodm_civ: fltarr(801,2), $  ; AODM error array for CIV
         ion: strarr(30), $
         wrest: dblarr(30), $   ;; Rest wavelength of transition
         wv_lim: dblarr(30,2), $  ;; Limits of integration in the spectrum
         flg_colm: intarr(30), $  ; 1 = Analysis, 2 = Lower limit, 4 = Upper limit
         instr: intarr(30), $     ; Instrument flag
         EW: fltarr(30), $        ; Rest EW (mA)
         sigEW: fltarr(30), $
         Ncolm: fltarr(30), $     ;AODM column density
         sigNcolm: fltarr(30), $
         b: fltarr(30), $       ;Doppler parameter (km/s)
         sigb: fltarr(30), $
         zabs: dblarr(30),  $      ; Redshift of line (tau-weighted)
         zsig: dblarr(30)  $       ; Error
         }

end
  
; Instrument 1=1bsic, 2=2asic, 4=2blif, 8=1alif, 16=1asic,
;     32=1blif, 64=2alif, 128+=STIS
         
