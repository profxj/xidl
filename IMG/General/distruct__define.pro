;+ 
; NAME:
; distruct__define
;    Version 1.1
;
; PURPOSE:
;    Defines the structure for direct imaging
;
; CALLING SEQUENCE:
;  tmp = {distruct}
;
; INPUTS:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   ??-2001 Written by JXP
;   11-June-2007 Modified by LKP for use with NIRC2
;
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro distruct__define

;  This routine defines the structure for direct images

  tmp = {distruct, $
         flg_anly: 0,$     ; Analysis flag 0=Don't Analyse
         frame: 0L,   $    ; FRAME Number
         naxis1: 0L,  $    ; NAXIS1 is helpful for determining sub-chip readouts, esp for NIRC2
         naxis2: 0L,  $    ; NAXIS2
         exp: 0.d,   $     ; Exposure time
         filter: ' ', $    ; Filter: U,B,V,R,I
         type: ' ',   $    ; ObjTyp: OBJ, STD, TWI, DRK, ZRO, DFT, SPS
         Obj: ' ', $       ; Object Name
         Equinox: 0.,$     ; EQUINOX
         Date: 0.0d,  $    ; Date of Obs
         UT: ' ',     $    ; UT
         RA: ' ',     $    ; RA
         DEC: ' ',    $    ; DEC
         TEL: ' ',    $    ; Telescope
         AM:   0.,   $     ; Airmass
         RMS:   0.,   $    ; RMS
         extinct:   0.,   $; Extinction (e.g. cirrus)
         CCD: ' ',    $    ; CCD
         namp: 0,     $    ; Namps
         coadds: 0L,  $    ; Number of coadds
         sampmode: 0L, $  ; Sampmode
         multisam: 0L, $   ; Multisam
         statsec: lonarr(4), $ ; Stat section
         gain: 0.,   $     ; Gain
         readno: 0., $     ; Read Noise
         satur: 0., $      ; Saturation
         rootpth: ' ',$    ; Path of the Root
         img_root: ' ',$   ; Root name (usually in Raw directory)
         med_raw: 0L, $    ; Median of the Raw image
         min_raw: 0L, $    ; Min of the Raw image
         max_raw: 0L, $    ; Max of the Raw image
         flg_ov:  0, $     ; OV FILE?  0=No, 1=Yes
         img_ov:  ' ', $   ; Name of OV file (with directory)
         med_ov:  0.,$     ; Median of OV file
         flg_drk: 0, $     ; Does necessary dark file exist? 0=NO, 1=Yes
         img_drk: ' ', $   ; Name of dark file which should be used with this image. (with directory)
         flg_drksub: 0, $  ; Has image been dark subtracted? 0=NO, 1=Yes
         img_drksub: ' ',$ ; Name of the dark subtracted (and coadd divided) file (with directory)
         med_drksub: 0.,$  ; Median of dark subtracted, coadded divided file
         flg_msk:  0,$     ; Mask FILE?  0=No, 1=Yes  
         img_msk:  ' ',$   ; Name of Mask file
         flg_skymsk:   0,$ ; Sky Mask FILE?  0=No, 1=Yes  
         img_skymsk:  ' ',$; Name of Sky Mask file
         flg_final: 0, $   ; Final File? 0=No, 1=OV, 2=Flat, 4=Gain
         img_final: ' '$   ; Name of Final img
         }

end
  
         
