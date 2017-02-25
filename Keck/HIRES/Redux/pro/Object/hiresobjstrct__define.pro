;+ 
; NAME:
; hiresobjstrct__define 
;   Version 1.1
;
; PURPOSE:
;    Defines the object structure which handles tracing and
;    extraction.
;
; CALLING SEQUENCE:
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
;   19-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------


pro hiresobjstrct__define

;  This routine defines the structure for individual object spectra

  tmp = {hiresobjstrct, $
         field: ' ', $
         order: 0L, $            ; Physical order #
         obj_id: ' ',        $   ; ID value (a=primary, b-z=serendip, x=NG)
         flg_anly: 0,      $     ; 0=No analy, 1=Traced, 2=Extracted, 4=Fluxed 
         exp: 0.d, $             ; Exposure time
         xcen: 0L, $             ; Column where obj was id'd
         ycen: 0., $
         skyrms: 0., $        
         spatial_fwhm: 0., $
         trace: fltarr(6000), $   ; Object trace
         fin_trc: fltarr(6000), $ ; Modified trace based on obj profile
         flg_aper: 0, $           ; 0=boxcar
         aper: fltarr(2), $       ; Width of aperture for masking (pixels)
         nrow: 0L, $
         box_wv: dblarr(6000), $  ; Box car extraction wavlengths
         box_fx: fltarr(6000), $  ; Box car extraction flux (electrons)
         box_var: fltarr(6000), $ ; Box car extraction variance (electrons)
         flg_sky: 0, $
         sky: fltarr(6000), $
         sky_wv: dblarr(6000), $
         skyshift: 0.d, $
         flg_optimal: 0, $
         npix: 0L, $
         opt_sigma: 0., $
         wave: dblarr(6000), $   ; Wavelengths for optimal extraction
         fx: fltarr(6000), $     ; Optimal fx
         var: fltarr(6000), $    ; <=0 :: rejected pix
         novar: fltarr(6000), $   ; <=0 :: rejected pix
         flg_flux: 0, $          ; 1=fnu
         flux: fltarr(6000), $   ; Fluxed data
         sig: fltarr(6000), $    ; Error in fluxed data
         nosig: fltarr(6000), $  ; Error in fluxed data not including obj counts
         date: 0.0d, $
         img_fil: ' ', $
         arc_fil: ' ', $
         UT: ' ' $
         }

end
  
         
