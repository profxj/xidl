pro fit2dstrct__define

;  This routine defines the Fit structure for a 2D surface fit

  tmp = {fit2dstrct, $
         func: ' ',  $       ; Name
         nx: 0,   $      ; Order number
         ny: 0,   $      ; Order number
         nrm: dblarr(2,2),   $      ; Normalization for x dimension
         lsig: 0., $
         hsig: 0., $
         niter: 0L, $
         minpt: 0L, $
         maxrej: 0L, $
         flg_rej: 0, $
         rms: 0., $
         ffit: ptr_new(/allocate_heap) $  ; 
         }

end
  
         
