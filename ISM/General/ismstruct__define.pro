;+ 
; NAME:
; ismstruct__define
;    Version 1.1
;
; PURPOSE:
;    Defines the structure for ISM (includes D/H)
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
pro ismstruct__define

;  This routine defines the structure for direct images

  tmpe = {elmstruct}
  tmpi = {ionstruct}

  tmp = {ismstruct, $
         ismfil: '',$
         target: '',$
         targ_ra: '',$
         targ_dec: '',$
         targ_vel: 0.0,$
         flg_mag: 0,$
         targ_mag: 0.0,$
         zabs: 0.0d,$
         NHI: 0.0d,$
         sigNHI: fltarr(2),$
         abndfil: '',$
         tab_fil: '', $
         flgD: 0,$
         DH: 0.0,$
         sigDH: 0.0,$
         flgFe: 0,$
         FeH: 0.0,$
         sigFeH: 0.0,$
         flgZn: 0,$
         ZnH: 0.0,$
         sigZnH: 0.0,$
         flgAlpha: 0,$
         Alpha: 0.0,$
         sigAlpha: 0.0,$
         flglw: 0,$
         lwfil: '',$
         lwwav: 0.0d,$
         lwvmn: 0.0,$
         lwvmx: 0.0,$
         lwfvel: 0.0,$
         lwfmm: 0.0,$
         lwfedg: 0.0,$
         lwftpk: 0.0,$
         flgCII: 0,$
         CII: 0.0,$
         sigCII: 0.0,$
         flgciv: 0.0,$
         civfil: '',$ 
         civwav: 0.0d,$
         civvmn: 0.0,$
         civvmx: 0.0,$
         civfvel: 0.0,$
         civfmm: 0.0,$
         civfedg: 0.0,$
         civftpk: 0.0,$
         civlwfdv: 0.0,$
         civlwfrto: 0.0,$
         civlwfnmm: 0.0,$
         civlwftvm: 0.0,$
         ebv: 0.0,$
         ffilt: 0,$
         fslit: 0,$
         srvy: 0, $
         srvy_mag: 0., $
         ref: '', $
         flgmtl: 0, $
         mtl: 0., $
         sigmtl: 0., $
         ndfil: 0,$
         Hfil: '',$
         Efil: '',$
         Ufil: '',$
         MBfil: '',$
         MRfil: '',$
         Xfil: '',$
         elm: replicate(tmpe,101),$  ; D is in 99
         XH: replicate(tmpe,101), $
         ion: replicate(tmpi,101) $
        }

end
  
         
