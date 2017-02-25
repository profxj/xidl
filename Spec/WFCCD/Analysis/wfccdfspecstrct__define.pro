pro wfccdfspecstrct__define

;  This routine defines the structure for individual object spectra

;  ??? mrb changed to not conflict with new idlspec2d conventions
  tmp1 = create_struct( $
    name = 'WFCCDZANS', $
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
    'vdisp_err'  , 0.0, $
    'vdispz'     , 0.0, $
    'vdispz_err' , 0.0, $
    'vdispchi2'  , 0.0, $
    'vdispnpix'  , 0.0, $
    'vdispdof'   ,  0L  $
   )

  tmp = {wfccdfspecstrct, $
         field: ' ', $
         slit_id: 0L, $          
         obj_id: ' ',        $   ; ID value (a=primary, b-z=serendip, x=NG)
         flg_anly: 0,      $     ;  0=No analy, 1=Extracted, 2=Fluxed 
         zans: tmp1, $           ; zans structure from SDSS package
         obj_type: 0, $
         xcen: 0L,$
         ycen: 0.,$
         trace: fltarr(5000),$
         xyimg: fltarr(2), $
         mag: 0., $              ; Usually R mag
         spec2d_fil: strarr(5), $
         phot_fil: ' ', $
         img_fil: ' ', $         ; Img file
         nexp: 0L, $
         wvmnx: fltarr(100,2), $ ; Wave min/max for each exposure
         texp: dblarr(100), $    ; t_exp for each exposure
         obj_fil: strarr(100), $
         npix: 0L, $
         flg_flux: 0, $          ; 1=fnu, 2=flambda
         wave: fltarr(4000), $
         fx: fltarr(4000), $
         var: dblarr(4000), $     ; <=0 :: rejected pix
         sdssname: ' ',$         ; INFO FROM SDSS FILES
         ra: 0.0d,$     
         dec: 0.0d,$
         rmag: 0.,$
         gr: 0.,$
         psfmodel: 0., $
         wfccd_pri:0.,$
         zsdss:0.0,$            ; SDSS REDSHIFT, =0 if no info
         berlind:0L, $           ; =1 if in Berlind group sample
         vgroup: 0.0$      ; group redshift
         }

end
  
         
