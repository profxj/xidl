pro mikestrct__define

;  This routine defines the structure for direct images
;  [v1.1]

  tmp = {mikestrct, $
         frame: 0,   $          ; FRAME Number
         flg_anly: 0,$          ; Analysis flag 0=Don't Analyse, 2=bias sub, 4=scatt light
         side: 0L, $            ; 1=blue; 2=red
         obj_id: 0L, $          ; Obj ID
         Obj: ' ', $            ; Object Name
         type: ' ',   $         ; ObjTyp: OBJ,STD,DRK,ZRO,IFLT,DFLT,ARC,MSK,IMG,TWI, MFLT, TFLT
         slit: 0., $            ; Slit
         setup: 0, $            ; Setup value
         exp: 0.d,   $          ; Exposure time
         colbin: 0L, $          ; Binning in column
         rowbin: 0L, $          ; Binning in row
         AM:   0.,   $          ; Airmass
         CCD: ' ',    $         ; CCD
         TEL: ' ',    $         ; Telescope
         ccdspeed: '', $        ; CCD Speed
         arclamp: 0, $          ; Arc Lamps: Binary  1=CuAr, 2=Xe, 4=HgNe
         arc_xyoff: fltarr(2),$ ; Offsets for Arc IMG due to thermal expansion
         qtzlamp: 0, $          ; QTZ Lamp: 0=Off, 1=On
         qtzmode: 0, $          ; 1=Normal 2=Milky
         gain: 0.,   $          ; Gain
         readno: 0., $          ; Read Noise
         date: 0.0d,  $         ; Date of Obs
         UT: 0L,     $          ; UT (s) [keyword UT-TIME]
         RA: ' ',     $         ; RA
         DEC: ' ',    $         ; DEC
         Equinox: 0.,$          ; EQUINOX
         rootpth: ' ',$         ; Path of the Root
         img_root: ' ',$        ; Root name (usually in Raw directory)
         flg_ov:  0, $          ; OV FILE?  0=No, 1=Yes  
         img_ov:  ' ', $        ; Name of OV file (with directory)
         flg_final: 0, $        ; Final File? 0=No
         img_final: ' ',$       ; Name of Final img
         ystrt: 0L, $           ; Column for initiating the trace
         arc_fil: ' ', $        ; Name of the Arc image file (fits)
         arc_img: ' ', $        ; Name of the Final Arc image file (fits)
         map_fil: ' ', $        ; Name of the Map image file (fits)
         flat_fil: ' ', $       ; Name of the Flat image file (fits)
         obj_fil: ' ' $         ; Object structure
         }

end
  
         
