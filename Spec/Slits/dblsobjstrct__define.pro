pro dblsobjstrct__define

;  This routine defines the structure for individual object spectra

  tmp = {dblsobjstrct, $
         field: ' ', $
         slit_id: 0L, $            ; Or order ID
         obj_id: ' ',        $      ; ID value (a=primary, b-z=serendip, x=NG)
         flg_anly: 0,      $  ;  0=No analy, 1=Extracted, 2=Fluxed 
         exp: 0.d, $
         xcen: 0L, $                ; Column where obj was id'd
         ycen: 0., $
         flg_aper: 0, $        ; 0=boxcar
         aper: fltarr(2), $            ; Widths of aperture (pixels)
         skyrms: 0., $
         trace: fltarr(5000), $
         npix: 0L, $
         wave: dblarr(5000), $
         fx: fltarr(5000), $
         var: fltarr(5000), $   ; <=0 :: rejected pix
         flg_flux: 0, $          ; 1=fnu
         flux: fltarr(5000), $   ; Fluxed data
         sig: fltarr(5000), $    ; Error in fluxed data
         date: 0.0d, $
         UT: ' ', $
         spec2d_fil: ' ', $
         img_fil: ' ', $
         img_xy: fltarr(2), $   ; xy pixel on the image
         slit_fil: ' ', $
         instr_strct: ' '$   ; e.g. wfccdstr fits file
         }

end
  
         
