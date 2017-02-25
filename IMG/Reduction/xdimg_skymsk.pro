;+ 
; NAME:
; xdimg_skymsk   
;   Version 1.1
;
; PURPOSE:
;    Builds a sky mask for each input image using xmkmask, an 
;      IDL gui built for this procedure.
;
; CALLING SEQUENCE:
;   
;   xdimg_skymask, struct, objimg, /ERASE, /DOOV, XSIZE=, YSIZE=
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;   objimg   -- integer array defining the images to mask
;
; RETURNS:
;
; OUTPUTS:
;   skymasks --  Builds sky masks and puts them into 'Masks/Sky/'
;
; OPTIONAL KEYWORDS:
;   erase - Erase pre-existing sky masks first
;   doov - Make OV files as necessary
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_skymsk, nght1_strct, objimg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-July-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_skymsk, struct, objimg, ERASE=erase, DOOV=doov, XSIZE=xsize, YSIZE=ysize

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'xdimg_skymsk, struct, objimg, /ERASE, /DOOV, XSIZE=, YSIZE= (v1.1)'
      return
  endif 
  
;  Optional Keywords

;  Overscan

  if keyword_set( DOOV ) then begin
      ovflt = where(struct[skyflt].flg_ov EQ 0, nov)
      if nov NE 0 then begin
          ovimg = struct[skyflt[ovflt]].rootpth + struct[skyflt[ovflt]].img_root
          xoverscan, ovimg, struct[skyflt[0]].ccd, ORDR=4
      endif
  endif


;  Determine images to mask

  if not keyword_set( ERASE ) then begin
      nomsk = where(struct[objimg].flg_skymsk EQ 0, Nnomsk)
      if Nnomsk EQ 0 then return else begin
          ; Mask
          xmkmask, struct[objimg[nomsk]].img_ov, /SKYMSK, XSIZE=xsize, YSIZE=ysize
          struct[objimg[nomsk]].flg_skymsk = 1
          struct[objimg[nomsk]].img_skymsk = 'Masks/Sky/sm_'+ $
            struct[objimg[nomsk]].img_root+'.gz'
      endelse
  endif else begin
      ; Mask
      xmkmask, struct[objimg].img_ov, /SKYMSK, XSIZE=xsize, YSIZE=ysize
      struct[objimg].flg_skymsk = 1
      struct[objimg].img_skymsk = 'Masks/Sky/sm_'+ $
        struct[objimg].img_root+'.gz'
  endelse

  ; Allow for a premature exit from xmkmask!
    
  for q=0,n_elements(objimg)-1 do begin
      if struct[objimg[q]].flg_skymsk EQ 1 then begin
          ; Does the file exist?
          a = findfile(struct[objimg[q]].img_skymsk, count=cnt)
          if cnt EQ 0 then struct[q].flg_skymsk = 0
      endif
  endfor

end
