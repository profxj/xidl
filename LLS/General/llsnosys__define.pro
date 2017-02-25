;+ 
; NAME:
; llsstruct__define
;    Version 1.1
;
; PURPOSE:
;    Defines the full structure for LLS
;
; CALLING SEQUENCE:
;  tmp = {ismstruct}
;
; INPUTS:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   7-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro llsnosys__define

  tmp = {llsnosys, $
         llsfil: '',$
         qso: '',$
         qso_ra: '',$
         qso_dec: '',$
         qso_zem: 0.0,$
         flg_mag: 0,$
         qso_mag: 0.0,$
         cldyfil: '', $
         flg_DH: 0, $
         DH: 0., $
         Nsys: 0, $
         NHI: 0., $
         sigNHI: fltarr(2), $
         NH: 0., $
         NHsig: fltarr(2), $
         zabs: 0., $
         flg_MH: 0., $
         MHave: 0., $
         MHsig: fltarr(2), $  ;; Used to be a single float
         vmn: 0.d, $
         vmx: 0.d, $
         fdelv: 0., $
         fmm: 0., $
         fedg: 0., $
         ftpk: 0., $
         tab_fil: '', $
         flglw: 0,$
         lwfil: '',$
         lwwav: 0.d,$
         lwvmn: 0.0,$
         lwvmx: 0.0,$
         lwfvel: 0.0,$
         lwfmm: 0.0,$
         lwfedg: 0.0,$
         lwftpk: 0.0,$
         srvy: 0, $
         srvy_mag: 0., $
         sdss_mjd: 0L, $
         sdss_plate: 0L, $
         sdss_fibid: 0L, $
         ref: '', $
         flgmtl: 0, $
         mtl: 0., $
         sigmtl: 0., $
         ndfil: 0,$
         Hfil: '',$
         Efil: '',$
         Ufil: '',$
         Xfil: '',$
         MBfil: '',$
         MRfil: '' $
;         elm: replicate(tmpe,101),$  ; D is in 99
;         XH: replicate(tmpe,101), $
;         ion:replicate(tmpi,101) $
        }

end
  
         
