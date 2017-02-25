;+ 
; NAME:
; dblsobjstrct__define
;   Version 1.1
;
; PURPOSE:
;  This routine creates a structure to describe the spectrum for an
;  object in a given slit
;
; CALLING SEQUENCE:
;   tmp = {dlbsobjstrct}
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
;  Written by JXP
;-
;----------------------------------------------------------------------
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
         gauss_sig: 0.0d,  $ ; gaussian sigma for profile fit
         skyrms: 0., $
         trace: fltarr(5000), $
         npix: 0L, $
         wave: dblarr(5000), $
         fx: fltarr(5000), $
         var: fltarr(5000), $   ; <=0 :: rejected pix
         novar: fltarr(5000), $   ; <=0 :: rejected pix
         sky:fltarr(5000), $
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
  
         
