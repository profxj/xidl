pro zoomstrct__define

;  This routine defines the structure for direct images

  tmp = {zoomstrct, $
         reg: lonarr(4),     $    ; Maps image onto entire window  
         level: 0L, $
         factor: 0., $
         flg: 0, $
         offset: lonarr(2), $
         centerpix: fltarr(2) $
         }

end
  
         
