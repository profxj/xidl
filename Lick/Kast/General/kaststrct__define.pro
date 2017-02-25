;+ 
; NAME:
; kaststrct__define   
;     Version 1.1
;
; PURPOSE:
;   Kast IDL structure
;
; CALLING SEQUENCE:
;  tmp = {kaststrct}
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro kaststrct__define

;  This routine defines the structure for Kast images
;  [v1.1]

  tmp = {kaststrct, $
         frame: 0,   $          ; FRAME Number
         flg_anly: 0,$          ; Analysis flag 0=Don't Analyse, 2=bias sub, 4=scatt light
         side: 0L, $            ; 1=blue; 2=red
         obj_id: 0L, $          ; Obj ID
         Obj: ' ', $            ; Object Name
         type: ' ',   $         ; ObjTyp: OBJ,STD,DRK,ZRO,QTZ,ARC,TWI,SLT,ACQ
         mode: 0, $             ; 0=Image, 1=Spec
         slit: 0., $            ; Slit
         splitter: ' ', $       ; Name
         setup: 0, $            ; Setup
         grising: ' ', $        ; Name of grism/grating
         exp: 0.d,   $          ; Exposure time
         AM:   0.,   $          ; Airmass
         CCD: ' ',    $         ; CCD
         PA: 0.,      $         ; PA
         TEL: ' ',    $         ; Telescope
         ccdspeed: '', $        ; CCD Speed
         arclamp: 0, $          ; Arc Lamps: Binary  1=Blue; 2=Red
         qtzlamp: 0, $          ; QTZ Lamp: 0=Off, 1=Blue; 2=SuperBlue
         gain: 0.,   $          ; Gain
         readno: 0., $          ; Read Noise
         date: 0.0d,  $         ; Date of Obs
         UT: ' ',     $         ; UT
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
         map_fil: ' ', $        ; Name of the Map image file (fits)
         flat_fil: ' ', $       ; Name of the Flat image file (fits)
         obj_fil: ' ' $         ; Object structure
         }

end
  
         
