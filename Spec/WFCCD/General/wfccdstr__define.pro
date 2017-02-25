pro wfccdstr__define

;  This routine defines the structure for direct images
;  [v1.1]

  tmp = {wfccdstr, $
         frame: 0,   $     ; FRAME Number
         flg_anly: 0,$     ; Analysis flag 0=Don't Analyse
         Obj: ' ', $       ; Object Name
         type: ' ',   $    ; ObjTyp: OBJ, STD, DRK, ZRO, FLT, ARC, MSK, ALG
         masknm: ' ', $    ; Mask Name: 'M1', 'M2', etc. or 'LONG'
         mask_id: 0L, $    ; Mask ID
         exp: 0.d,   $     ; Exposure time
         filter: ' ', $    ; Filter: U,B,V,R,I, C
         filtpos: 0, $     ; Filter Position
         grism: ' ', $     ; Grism: BG, RG, NO
         aperpos: 0, $     ; Aperture position
         casspos: 0., $    ; CASS pos angle
         AM:   0.,   $     ; Airmass
         CCD: ' ',    $    ; CCD
         TEL: ' ',    $    ; Telescope
         gain: 0.,   $     ; Gain
         readno: 0., $     ; Read Noise
         date: 0.0d,  $    ; Date of Obs
         UT: ' ',     $    ; UT
         RA: ' ',     $    ; RA
         DEC: ' ',    $    ; DEC
         Equinox: 0.,$     ; EQUINOX
         rootpth: ' ',$    ; Path of the Root
         img_root: ' ',$   ; Root name (usually in Raw directory)
         flg_ov:  0, $     ; OV FILE?  0=No, 1=Yes  
         img_ov:  ' ', $   ; Name of OV file (with directory)
         flg_msk:  0,$     ; Mask FILE?  0=No, 1=Yes  
         img_msk:  ' ',$   ; Name of Mask file
         flg_final: 0, $   ; Final File? 0=No
         img_final: ' ',$  ; Name of Final img
         nslits: 0, $      ; Num of slits
         ystrt: 0L, $      ; Column for initiating the trace
         msk_fil: ' ', $   ; Name of the Mask info file (fits)  'Masks/q0026..'
         slit_fil: ' ', $  ; Name of the Slit info file (fits)
         arc_fil: ' ', $   ; Name of the Arc image file (fits)
         map_fil: ' ', $   ; Name of the Map image file (fits)
         flat_fil: ' ', $  ; Name of the Flat image file (fits)
         obj_fil: ' ', $   ; Object structure
         rotated:0L, $     ; 1L if we are setting rotation by hand
         theta:0., $       ;   angle of rotation
         shift:0. $        ;   shift of image
         }

end
  
         
