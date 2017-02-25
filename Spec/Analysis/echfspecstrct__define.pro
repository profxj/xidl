;+ 
; NAME:
; echfspecstrct__define
;   Version 1.1
;
; PURPOSE:
;  Creates a structure for echelle spectroscopy that will hold the
;  wavelength, flux and error arrays.  Also includes the ZANS
;  structure which is useful for using SDSS redshift identification.
; 
;
; CALLING SEQUENCE:
;   tmp = {echfspecstrct}
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro echfspecstrct__define

;  This routine defines the structure for individual object spectra

  tmp1 = create_struct( $
    name = 'ZANS', $
    'class'      ,  '', $
    'subclass'   ,  '', $
    'z'          , 0.0, $
    'z_err'      , 0.0, $
    'rchi2'      , 0.0, $
    'dof'        ,  0L, $
    'rchi2diff'  , 0.0, $
    'tfile'      ,  '', $
    'tcolumn'    , lonarr(20) - 1L, $
    'npoly'      ,  0L, $
    'theta'      , fltarr(20), $
    'vdisp'      , 0.0, $
    'vdisp_err'  , 0.0  $
   )

  tmp = {echfspecstrct, $
         field: ' ', $
         obj_id: ' ',        $   ; ID value (a=primary, b-z=serendip, x=NG)
         flg_anly: 0,      $     ;  0=No analy, 1=Extracted, 2=Fluxed 
         zans: tmp1, $           ; zans structure from SDSS package
         obj_type: 0, $
         xyimg: fltarr(2), $     ; xy pix of original image
         mag: 0., $              ; Usually R mag
         phot_fil: ' ', $
         img_fil: ' ', $         ; Img file
         nexp: 0L, $
         wvmnx: dblarr(100,2), $ ; Wave min/max for each exposure
         texp: dblarr(100), $    ; t_exp for each exposure
         obj_fil: strarr(100), $
         flg_fin: 0, $           ; 
         flg_flux: 0, $          ; 1=fnu, 2=flambda
         npix: lonarr(50), $     ; ORDERS
         wave: dblarr(5000,50), $
         fx: fltarr(5000,50), $
         var: dblarr(5000,50), $     ; <=0 :: rejected pix
         novar:dblarr(5000, 50), $
         sky: dblarr(5000, 50) $
         }

end
  
         
