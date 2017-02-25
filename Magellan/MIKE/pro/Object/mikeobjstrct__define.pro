pro mikeobjstrct__define

;  This routine defines the structure for individual object spectra

  tmp = {mikeobjstrct, $
         field: ' ', $
         order: 0L, $            ; Physical order #
         obj_id: ' ',        $   ; ID value (a=primary, b-z=serendip, x=NG)
         flg_anly: 0,      $     ; 0=No analy, 1=Traced, 2=Extracted, 4=Fluxed 
         exp: 0.d, $             ; Exposure time
         xcen: 0L, $             ; Column where obj was id'd
         ycen: 0., $
         skyrms: 0., $        
         spatial_fwhm: 0., $
         trace: fltarr(5000), $   ; Object trace
         fin_trc: fltarr(5000), $ ; Modified trace based on obj profile
         flg_aper: 0, $           ; 0=boxcar
         aper: fltarr(2), $       ; Width of aperture for masking (pixels)
         nrow: 0L, $
         box_fwhm: fltarr(5), $
         box_wv: dblarr(5000), $  ; Box car extraction wavlengths
         box_fx: fltarr(5000), $  ; Box car extraction flux (electrons)
         box_var: fltarr(5000), $ ; Box car extraction variance (electrons)
         flg_sky: 0, $
         sky: fltarr(5000), $
         sky_wv: dblarr(5000), $
         skyshift: 0.d, $
         flg_optimal: 0, $
         npix: 0L, $
         wave: dblarr(5000), $   ; Wavelengths for optimal extraction
         fx: fltarr(5000), $     ; Optimal fx
         var: fltarr(5000), $    ; <=0 :: rejected pix
         novar: fltarr(5000), $  ; added by jm08apr15nyu
         flg_flux: 0, $          ; 1=fnu
         flux: fltarr(5000), $   ; Fluxed data
         sig: fltarr(5000), $    ; Error in fluxed data
         nosig: fltarr(5000), $  ; added by jm08apr15nyu
         date: 0.0d, $
         img_fil: ' ', $
         arc_fil: ' ', $
         UT: ' ' $
         }

end
  
         
