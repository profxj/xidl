;+ 
; NAME:
; esi_lwdcombimg   
;     Version 1.0
;
; PURPOSE:
;    Process a data frame
;
; CALLING SEQUENCE:
;   
;  esi_lwdcombimg, esi, indx, /DFLAT, /REDDOV
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
;   esi_lwdcombimg, esi, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function esi_lwdcombimg, esi, indx, VAR=var, SILENT=silent, $
                         IMGINDX=imgindx, SKY=sky

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'img = esi_lwdcombimg( esi, indx, VAR= )[v1.0]'
      return, -1
  endif 
  
;  Optional Keywords

; NIMG
  nimg = n_elements(indx)

; NIMG = 1
  if nimg LT 2 then begin
      if not keyword_set( SILENT ) then $
        print, 'esi_lwdcombimg: Only 1 image, returning it'
      img = xmrdfits(esi[indx].img_final, 0, /silent)
      if arg_present( VAR ) then var = xmrdfits(esi[indx].img_final, 1, /silent)
  endif

; NIMG = 2
  if NIMG EQ 2 then img = esi_addtwo(esi, indx, VAR=var, /SCALE, SKY=SKY)

; NIMG > 2
  if NIMG GT 2 then begin
      ;; IMG
      xcombine, esi[indx].img_final, img, FCOMB=2, SCALE=esi[indx].exp, $
        GAIN=esi[indx[0]].gain, RN=esi[indx[0]].readno,IMGINDX=imgindx
      ;; VAR
      if arg_present( VAR ) then $
        xcombine, esi[indx].img_final, var, FCOMB=2, SCALE=esi[indx].exp, $
        GAIN=esi[indx[0]].gain, RN=esi[indx[0]].readno,IMGINDX=1L
  endif
  
  return, img
end
