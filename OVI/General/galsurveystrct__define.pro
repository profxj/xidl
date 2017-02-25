;+ 
; NAME:
; galsurveystrct__define
;    Version 1.1
;
; PURPOSE:
;  This routine defines the structure for galaxies in the OVI survey
;
; CALLING SEQUENCE:
;   tmp = {galsurveystrct}
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
;   19-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro galsurveystrct__define

;  This routine defines the structure for individual object spectra

  tmp = {galsurveystrct, $
         field: ' ', $
         id: 0L, $
         obj_id: ' ',        $   ; ID value (a=primary, b-z=serendip, x=NG)
         flg_anly: 0,      $     ;  0=No analy, 1=Phot, 2=Spec, 4=redshift
         flg_survey: 0,      $     ;  0=Not incluced; 1=Included
         obj_type: 0, $
         mag: fltarr(10), $        ; Mag
         magerr: fltarr(10), $     ; Error in Mag
         filter: strarr(10), $     ; Filter
         img_fil: strarr(10), $    ; Img file
         xypix: fltarr(2), $
         ra: 0.d, $
         dec: 0.d, $
         area: 0., $               ; sqr''
         stargal: 0., $            ; Star/galaxy classifier
         gal_type: ' ', $
         gal_coeff: fltarr(10), $  ; 0-3 are the eigen coeff; 4 = Late value
         z: 0.d, $
         vcirc: 0., $
         fspec_fil: strarr(10) $   ; Spectrocscopy files
         }

end
  
