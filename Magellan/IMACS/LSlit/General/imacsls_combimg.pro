;+ 
; NAME:
; imacsls_combimg   
;     Version 1.1
;
; PURPOSE:
;    Process a data frame
;
; CALLING SEQUENCE:
;  img = imacsls_combimg(imacsls, indx, IVAR=, /SILENT, IMGINDX=, /SKY)
;   
; INPUTS:
;   imacsls -  IMACS structure
;   indx    -  Index values
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  IVAR=  -- Inverse variance image
;
; COMMENTS:
;
; EXAMPLES:
;   imacsls_combimg, imacsls, indx
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function imacsls_combimg, imacsls, indx, IVAR=ivar, SILENT=silent, $
                         IMGINDX=imgindx, SKY=sky

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'img = imacsls_combimg( imacsls, indx, IVAR=, /SILENT, IMGINDX= )[v1.1]'
      return, -1
  endif 
  
;  Optional Keywords
  if not keyword_set( IMGINDX ) then imgindx = 0L
  if not keyword_set( VARINDX ) then varindx = 1L

; NIMG
  nimg = n_elements(indx)

; NIMG = 1
  if nimg LT 2 then begin
      if not keyword_set( SILENT ) then $
        print, 'imacsls_combimg: Only 1 image, returning it'
      img = xmrdfits(imacsls[indx].img_final, IMGINDX, /silent)
      if arg_present( IVAR ) then $
        ivar = xmrdfits(imacsls[indx].img_final, 1, /silent)
  endif

; NIMG = 2
  if NIMG EQ 2 then img = imacsls_addtwo(imacsls, indx, IVAR=ivar, /SCALE, SKY=sky)

; NIMG > 2
  if NIMG GT 2 then begin
      ;; Add IMG with weighting
      x_addimg, imacsls[indx].img_final, img, ivar, $
        SCALE=mean(imacsls[indx].exp)/imacsls[indx].exp, $
        WEIGHT=imacsls[indx].exp, IMGINDX=imgindx, VARINDX=varindx
  endif
  
  return, img
end
