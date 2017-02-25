pro esistrct__define

;  This routine defines the structure for direct images
;  [v1.1]

  tmp = {esistrct, $
         frame: 0,   $          ; FRAME Number
         flg_anly: 0,$          ; Analysis flag 0=Don't Analyse, 2=bias sub, 4=scatt light
         obj_id: 0L, $          ; Obj ID
         Obj: ' ', $            ; Object Name
         type: ' ',   $         ; ObjTyp: OBJ,STD,DRK,ZRO,IFLT,DFLT,ARC,MSK,IMG,TWI
         slit: 0., $            ; Slit
         exp: 0.d,   $          ; Exposure time
         imfilt: ' ', $         ; Image Filter: U,B,V,R,I, C
         mode: 0L, $            ; Mode:  0=Image, 1=LowD, 2=Echellette
         cbin: 0L, $            ; Column bin
         rbin: 0L, $            ; Row bin
         AM:   0.,   $          ; Airmass
         CCD: ' ',    $         ; CCD
         TEL: ' ',    $         ; Telescope
         namp: 0, $             ; Number of Amplifiers
         arclamp: 0, $          ; Arc Lamps: Binary  1=CuAr, 2=Xe, 4=HgNe
         qtzlamp: 0, $          ; QTZ Lamp: 0=Off, 1=On
         rotmode: 0, $          ; Rotation mode: 0=Stationary, 
         ccdspeed: '', $        ; CCD Speed
         gain: 0.,   $          ; Gain
         readno: 0., $          ; Read Noise
         date: 0.0d,  $         ; Date of Obs
         UT: ' ',     $         ; UT
         RA: ' ',     $         ; RA
         DEC: ' ',    $         ; DEC
         Equinox: 0.,$          ; EQUINOX
         refordr: 0L, $         ; Reference order
         rootpth: ' ',$         ; Path of the Root
         img_root: ' ',$        ; Root name (usually in Raw directory)
         flg_ov:  0, $          ; OV FILE?  0=No, 1=Yes  
         img_ov:  ' ', $        ; Name of OV file (with directory)
         flg_final: 0, $        ; Final File? 0=No
         img_final: ' ',$       ; Name of Final img
         ystrt: 0L, $           ; Column for initiating the trace
         arc_fil: ' ', $        ; Name of the Arc image file (fits)
         map_fil: ' ', $        ; Name of the Map image file (fits)
         flat_fil: ' ', $       ; Name of the Flat image file (fits)
         twiflat_fil: ' ', $     ; Name of the Twilight Flat image file (fits)
         obj_fil: ' ' $         ; Object structure
         }

end
  
         
