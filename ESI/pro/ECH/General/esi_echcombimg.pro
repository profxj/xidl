;+ 
; NAME:
; esi_echcombimg   
;     Version 1.0
;
; PURPOSE:
;    Process a data frame
;
; CALLING SEQUENCE:
;   
;  esi_echcombimg, esi, indx, /DFLAT, /REDDOV
;
; INPUTS:
;   esi     -  ESI structure
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
;   esi_echcombimg, esi, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function esi_echcombimg, esi, indx, VAR=var, SILENT=silent, $
                         IMGINDX=imgindx

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'img = esi_echcombimg( esi, indx, VAR= )[v1.0]'
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
        print, 'esi_echcombimg: Only 1 image, returning it'
      img = xmrdfits(esi[indx].img_final, IMGINDX, /silent)
      if arg_present( VAR ) then $
        var = xmrdfits(esi[indx].img_final, 1, /silent)
  endif

; NIMG = 2
  if NIMG EQ 2 then img = esi_echaddtwo(esi, indx, VAR=var, /SCALE)

; NIMG > 2
  if NIMG GT 2 then begin
      ;; Add IMG with weighting
      x_addimg, esi[indx].img_final, img, var, $
        SCALE=mean(esi[indx].exp)/esi[indx].exp, $
        WEIGHT=esi[indx].exp, IMGINDX=imgindx, VARINDX=varindx
  endif
  
  return, img
end
