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
pro llsstruct__define

;  This routine defines the structure for direct images

  tmpe = {elmstruct}
  tmpi = {ionstruct}
  tmpxh = {llsxhstruct}

  tmps = { sysstrct, $
           name: '', $
           zabs: 0.d, $
           NHI: 0., $
           NHIsig: fltarr(2), $
           NH: 0., $
           NHsig: fltarr(2), $
           logx: 0., $  ;; actually log(1-x)
           sig_logx: fltarr(2), $ ;; minus, plus
           b: 0., $
           bsig: 0., $
           U: 0., $
           Usig: fltarr(2), $
           abndfil: '', $
           flg_alpha: 0, $
           alphaH: 0., $
           ;sig_alphaH: 0., $
           sig_alphaH: fltarr(2), $
           flgFe: 0, $
           FeH: 0., $
           ;sig_FeH: 0., $
           sig_FeH: fltarr(2), $
           ndfil: 0,$
           Hfil: '',$
           Efil: '',$
           Ufil: '',$
           Xfil: '',$
           MBfil: '',$
           MRfil: '',$
           vpfil: '',$
           tab_fil: '', $
           elm: replicate(tmpe,101),$ ; D is in 99
	   XH: replicate(tmpxh,101), $ ; assymetric errors
           ion: replicate(tmpi,101), $
           flg_low: 0}

  tmp = {llsstruct, $
         llsfil: '',$
         qso: '',$
         qso_ra: '',$
         qso_dec: '',$
         qso_zem: 0.0,$
         flg_mag: 0,$
         qso_mag: 0.0,$
         systems: replicate(tmps,10), $
         cldyfil: '', $
         flg_DH: 0, $
         DH: 0., $
         Nsys: 0, $
         NHI: 0., $
         sigNHI: fltarr(2), $
         NH: 0., $
         NHsig: fltarr(2), $
         zabs: 0., $
         flg_MH: 0, $
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
         MRfil: '',$
         ion:replicate(tmpi,101) $
        }

end
  
         
