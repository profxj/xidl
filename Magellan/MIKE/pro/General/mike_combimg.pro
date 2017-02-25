;+ 
; NAME:
; mike_combimg   
;     Version 1.0
;
; PURPOSE:
;    Process a data frame
;
; CALLING SEQUENCE:
;   
;  mike_combimg, mike, indx, /DFLAT, /REDDOV
;
; INPUTS:
;   mike     -  ESI structure
;   indx    -  Index values
;
; RETURNS:
;
; OUTPUTS:
;  Fully combimgessed image
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_combimg, mike, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mike_combimg, mike, indx, IVAR=ivar, SILENT=silent, $
                         IMGINDX=imgindx, SKY=sky

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'img = mike_combimg( mike, indx, IVAR= )[v1.0]'
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
        print, 'mike_combimg: Only 1 image, returning it'
      img = xmrdfits(mike[indx].img_final, IMGINDX, /silent)
      if arg_present( IVAR ) then $
        ivar = xmrdfits(mike[indx].img_final, 1, /silent)
  endif

; NIMG = 2
  if NIMG EQ 2 then img = mike_addtwo(mike, indx, IVAR=ivar, /SCALE, SKY=sky)

; NIMG > 2
  if NIMG GT 2 then begin
      ;; Add IMG with weighting
      x_addimg, mike[indx].img_final, img, ivar, $
        SCALE=mean(mike[indx].exp)/mike[indx].exp, $
        WEIGHT=mike[indx].exp, IMGINDX=imgindx, VARINDX=varindx
  endif
  
  return, img
end
