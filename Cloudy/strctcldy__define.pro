pro strctcldy__define

;  This routine defines the structure for the CLOUDY grid

  tmp = {strctcldy, $
         z: 0.d,$          ; Redshift
         NHI: 0.d,   $     ; N(HI)
         FeH: 0.d, $       ; Metallicity
         U: 0.d,   $       ; Ionization Parameter
         T: 0. ,   $       ; Temperature
         nH:   0.d,   $    ; volume density
         Jnu: 0.d,    $    ; Jnu
         Spec: ' ',    $   ; Spectrum shape (HM = Haardt-Madau)
         flg: 0, $         ; Flag (0 = No output, 1 = Output)
         X: fltarr(31,31) $ ; Ionic ratios 
         }

end
  
         
