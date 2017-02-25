;+ 
; NAME:
; xdimg_objmsk   
;   Version 1.0
;
; PURPOSE:
;    Builds a sky mask for each input image using xmkmask, an 
;      IDL gui built for this procedure.
;
; CALLING SEQUENCE:
;   xdimg_skymask, struct, objimg, XSIZE=, YSIZE=, IMSK=
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
;   XSIZE, YSIZE -- Size of GUI
;
; OPTIONAL OUTPUTS:
;  IMSK -- Image of mask
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_objmsk, nght1_strct, objimg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_objmsk, struct, objimg, XSIZE=xsize, YSIZE=ysize, IMSK=imsk

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'xdimg_objmsk, struct, objimg, /ERASE, /DOOV (v1.0)'
      return
  endif 
  
;  Optional Keywords

;  Determine images to mask

      ; Mask
  xmkmask, struct[objimg].img_final, XSIZE=xsize, YSIZE=ysize, $
    IMSK=imsk
  struct[objimg].flg_msk = 1
  struct[objimg].img_msk = 'Masks/mk_'+ $
    struct[objimg].img_root+'.gz'

  ; Allow for a premature exit from xmkmask!
    
;  for q=0,n_elements(objimg)-1 do begin
;      if struct[objimg[q]].flg_objmsk NE 0 then begin
;          ; Does the file exist?
;          a = findfile(struct[objimg[q]].img_objmsk, count=cnt)
;          if cnt EQ 0 then struct[q].flg_objmsk = 0
;      endif
;  endfor

end
