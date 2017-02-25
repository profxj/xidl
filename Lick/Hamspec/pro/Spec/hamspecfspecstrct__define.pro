;+
; NAME:
; hamspec_1dspec
;    Version 1.1
;
; PURPOSE:
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
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-May-2013 Slightly modified from hiresfspecstrct__define by ENK
;-
pro hamspecfspecstrct__define

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

  tmp = {hamspecfspecstrct, $
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
         npix: lonarr(160), $     ; ORDERS
         phys_ordr: lonarr(160), $     ; ORDERS
         wave: dblarr(6000L,160), $
         fx: fltarr(6000L,160), $
         var: dblarr(6000L,160), $     ; <=0 :: rejected pix
         novar: dblarr(6000L,160) $
         }

end
  
         
