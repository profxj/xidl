pro mslitstrct__define

; Version 1.0

  tmp = {mslitstrct, $
         flg_anly: 0, $
         id: 0L,        $      ; ID number
         field: ' ',   $      ; Field name
         length: 0., $
         width: 0., $
         pa: 0., $
         arcpix: 0., $
         xyqso: fltarr(2),    $       ; x offset from the QSO (arcsec)
         xpos: 0., $ ; x position of silt
         nobj: 0, $  ; Number of objects
         ypos: fltarr(10), $  ; y pos of obj relative to slit edges (percentage)
         yedg: fltarr(2), $   ; Slit edges relative to slit mask cen (pix)
         yedg_flt: fltarr(2),$; y-Edges  of slit on CCD (Undistorted FLAT)
         yedg_orig: fltarr(5000,2), $  ; Array of slit edges in orig (1=top)
         yedg_sky: fltarr(5000,2)  $  ; Array of slit edges for sky sub
         }

end
  
         
