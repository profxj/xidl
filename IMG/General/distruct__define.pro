pro distruct__define

;  This routine defines the structure for direct images

  tmp = {distruct, $
         flg_anly: 0,$     ; Analysis flag 0=Don't Analyse
         frame: 0L,   $    ; FRAME Number
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
         statsec: lonarr(4), $ ; Stat section
         gain: 0.,   $     ; Gain
         readno: 0., $     ; Read Noise
         satur : 0., $     ; Saturation
         rootpth: ' ',$    ; Path of the Root
         img_root: ' ',$   ; Root name (usually in Raw directory)
         med_raw: 0L, $    ; Median of the Raw image
         min_raw: 0L, $    ; Min of the Raw image
         max_raw: 0L, $    ; Max of the Raw image
         flg_ov:  0, $     ; OV FILE?  0=No, 1=Yes  
         img_ov:  ' ', $   ; Name of OV file (with directory)
         med_ov:  0.,$     ; Median of OV file
         flg_msk:  0,$     ; Mask FILE?  0=No, 1=Yes  
         img_msk:  ' ',$   ; Name of Mask file
         flg_skymsk:   0,$ ; Sky Mask FILE?  0=No, 1=Yes  
         img_skymsk:  ' ',$; Name of Sky Mask file
         flg_final: 0, $   ; Final File? 0=No, 1=OV, 2=Flat, 4=Gain
         img_final: ' '$   ; Name of Final img
         }

end
  
         
