;+ 
; NAME:
; xdimg_updstr   
;   Version 1.1
;
; PURPOSE:
;    Looks around for OV and other images and resets the appropriate
;   flags.
;
; CALLING SEQUENCE:
;  xdimg_updstr, strct, /REDOMED
;
; INPUTS:
;   strct -- Direct Image strcure DISTRUCT
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
;   xdimg_updstr, struct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_updstr, struct, REDOMED=redomed

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'xdimg_updstr, struct, /REDOMED (v1.1)'
      return
  endif 
  
  nimg = n_elements(struct)

;  Loop on Indv images

  for q=0,nimg-1 do begin

      ;; OV
      ovfil = 'OV/ov_'+struct[q].img_root
      a = findfile(ovfil, count=cnt)
      if cnt NE 0 then begin
          struct[q].flg_ov = 1
          struct[q].img_ov = ovfil
          ; Median
          if struct[q].med_ov EQ 0 or keyword_set( REDOMED ) then begin
              xoverscan, ovfil, struct[q].ccd, medov, /MEDONLY
              struct[q].med_ov = medov
          endif
      endif else struct[q].flg_ov = 0
     
      ;; Mask
      mskfil = 'Masks/mk_'+struct[q].img_root+'.gz'
      a = findfile(mskfil, count=cnt)
      if cnt NE 0 then begin
          struct[q].flg_msk = 1
          struct[q].img_msk = mskfil
      endif else struct[q].flg_msk = 0

      ;; SkyMask
      skyfil = 'Masks/Sky/sm_'+struct[q].img_root+'.gz'
      a = findfile(skyfil, count=cnt)
      if cnt NE 0 then begin
          struct[q].flg_skymsk = 1
          struct[q].img_skymsk = skyfil
      endif else struct[q].flg_skymsk = 0
      
      ;; FINAL
      finfil = 'Final/f_'+struct[q].img_root
      a = findfile(finfil, count=cnt)
      if cnt NE 0 then begin
          struct[q].flg_final = 1
          struct[q].img_final = finfil
      endif else struct[q].flg_final = 0

  endfor
      
  if not keyword_set( NOFILE ) then write_dimgstr, struct

  print, 'All done updating!'

end
