pro extrctstrct__define

;  This routine defines the exraction structure

  tmp1 = { fitstrct }
  tmp = {extrctstrct, $
         ovr: fltarr(4), $ ; OV Region 
         aper: fltarr(2),  $  ; Aperture
         cline: 0L,   $      ; Line defining the aperture
         skyreg: PTR_NEW(/allocate_heap),   $  ; Sky Regions
         skyfstr: tmp1, $ ; Sky fit structure
         trace: PTR_NEW(/allocate_heap) $  ; f-value   
         }

end
  
         
