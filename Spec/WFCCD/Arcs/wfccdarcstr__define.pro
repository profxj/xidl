pro wfccdarcstr__define

;  This routine defines the structure for direct images

  tmp1 = { fitstrct }
  tmp = {wfccdarcstr, $
         flg_anly: 0, $            ; Good fit
         cent: 0., $               ; Region where arc spect is determined
         wave: dblarr(2048), $     ; Wavelength array
         spec: fltarr(2048), $     ; Arc spectrum
         fit: tmp1 $               ; Fit structure
         }

end
  
         
