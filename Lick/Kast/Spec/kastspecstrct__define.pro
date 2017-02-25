pro kastspecstrct__define

;  This routine defines the structure for individual object spectra

;    name = 'ZANS', $
  tmp1 = { xzans, $
           class:  '', $
           subclass: '', $
           z:          0.0, $
           z_err:      0.0, $
           rchi2:      0.0, $
           dof:        0L, $
           rchi2diff:  0.0, $
           tfile:      '', $
           tcolumn:    lonarr(20) - 1L, $
           npoly:      0L, $
           theta:      fltarr(20), $
           vdisp:      0.0, $
           vdisp_err:  0.0  $
         }

;    'rchi2diff'  , 0.0, $
;    'tfile'      ,  '', $
;    'tcolumn'    , lonarr(20) - 1L, $
;    'npoly'      ,  0L, $
;    'theta'      , fltarr(20), $
;    'vdisp'      , 0.0, $
;    'vdisp_err'  , 0.0  $
;   )

  tmp = {kastspecstrct, $
         field: ' ', $
         slit_id: 0L, $
         obj_id: ' ',        $   ; ID value (a=primary, b-z=serendip, x=NG)
         flg_anly: 0,      $     ;  0=No analy, 1=Extracted, 2=Fluxed 
         zans: tmp1, $           ; zans structure from SDSS package
         obj_type: 0, $
         xyimg: fltarr(2), $     ; xy pix of original image
         mag: 0., $              ; Usually R mag
         phot_fil: ' ', $
         img_fil: ' ', $         ; Img file
         nexp: 0L, $
         wvmnx: fltarr(100,2), $ ; Wave min/max for each exposure
         texp: dblarr(100), $    ; t_exp for each exposure
         obj_fil: strarr(100), $
         npix: 0L, $
         flg_flux: 0, $          ; 1=fnu, 2=flambda
	 flux: fltarr(2000), $
         wave: fltarr(2000), $
         fx: fltarr(2000), $
         var: dblarr(2000), $     ; <=0 :: rejected pix
         sig: fltarr(2000) $
         }

end
  
         
