;+ 
; NAME:
; dlastruct__define
;    Version 1.1
;
; PURPOSE:
;    Defines the structure for DLA
;
; CALLING SEQUENCE:
;  tmp = {dlastruct}
;
; INPUTS:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   1-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro dlastruct__define

;  This routine defines the structure for direct images

  tmpe = {elmstruct}
  tmpi = {ionstruct}

  tmp = {dlastruct, $
         dlafil: '',$
         qso: '',$
         qso_ra: '',$
         qso_dec: '',$
         qso_zem: 0.0,$
         flg_QSOmag: 0,$
         qso_mag: 0.0,$
         zabs: 0.0d,$
         NHI: 0.0d,$
         sigNHI: fltarr(2),$
         abndfil: '',$
         tab_fil: '', $
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
         qso_ebv: 0.0,$
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
         Xfil: '',$
         Ffil: '',$
         elm: replicate(tmpe,101),$
         XH: replicate(tmpe,101), $
         ion: replicate(tmpi,101), $
         vpfit_fil: '', $  ;; File for VPFIT.  Usually for CI and/or H2
         flg_H2: 0, $
         H2: 0., $  ;; Column of H_2
         sig_H2: 0., $
         flg_CI: 0, $  
         CI: 0., $  ;; N(C^0)
         sig_CI: 0., $
         sdss_mjd: 0L, $
         sdss_plate: 0L, $
         sdss_fibid: 0L $
        }

end
  
         
