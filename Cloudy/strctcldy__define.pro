;+ 
; NAME:
; strctcldy__define
;  V1.1
;
; PURPOSE:
;    Named IDL Structure for Cloudy output
;
; CALLING SEQUENCE:
;   tmp = {strctcldy}
;
; REVISION HISTORY:
;   31-May-2001 Written by JXP
;-
;------------------------------------------------------------------------------
pro strctcldy__define

;  This routine defines the structure for the CLOUDY grid

  tmp = {strctcldy, $
         z: 0.d,$          ; Redshift
         NHI: 0.d,   $     ; N(HI)
         FeH: 0.d, $       ; Metallicity
         U: 0.d,   $       ; Ionization Parameter
         Tval: 0. ,   $       ; Stopping Temperature
         Heat: 0., $       ; Heat source  (log erg s^-1)
         Tgas: fltarr(2) ,   $  ; Gas temperature [HI and HII]
         nH:   0.d,   $    ; volume density
         Jnu: 0.d,    $    ; Jnu
         Spec: ' ',    $   ; Spectrum shape (HM = Haardt-Madau)
         flg: 0, $         ; Flag (0 = No output, 1 = Output)
         X: fltarr(31,31) $ ; Ionic ratios 
         }

end
  
         
