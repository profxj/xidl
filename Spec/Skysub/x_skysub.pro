;+ 
; NAME:
; x_skysub   
;   Version 1.1
;
; PURPOSE:
;    Sky subtracts an input image
;
; CALLING SEQUENCE:
;   newimg = x_skysub(img, SKY=, REG=, CREG=, SKYFSTR=,
;                       TRACE=, /NOREJ)
;
; INPUTS:
;   img        - OV subtracted image
;
; RETURNS:
;   newimg     - Sky subtracted image
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;    REG       -  Sky region
;    CREG      -  Column used to define sky regions
;    TRACE     -  1-D trace to the object
;    SKYFSTR   -  Fit strcture for the sky
;    /NOREJ    -  Do not reject bad points
;
; OPTIONAL OUTPUTS:
;    sky       - 2-D image of the sky fit
;
; COMMENTS:
;
; EXAMPLES:
;   fimg = x_skysub(img, REG=skyreg)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Dec-2001 Written by JXP
;   07-Feb-2002 Revised to take fit structure
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_skysub, img, sky, REG=reg, CREG=creg, SKYFSTR=skyfstr, $
                   TRACE=trace, NOREJ=norej, SKYRMS=skyrms


;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'nwimg = x_skysub(img, [sky], REG=, CREG=, SKYFSTR=,' 
    print, '        TRACE=, NOREJ=, SKYRMS) [v1.1]'
    return, -1
  endif 

;  Optional Keywords

  if not keyword_set( SKYFSTR ) then begin
      skyfstr = { fitstrct }
      skyfstr.func = 'POLY'
      skyfstr.nord = 2
      if not keyword_set( NOREJ ) then begin
          skyfstr.hsig = 2.
          skyfstr.lsig = 2.
          skyfstr.niter = 3
      endif
  endif else begin
      if keyword_set( NOREJ ) then skyfstr.niter = 0
  endelse


  sz = size(img, /dimensions)
  if keyword_set( REG ) AND not keyword_set( CREG ) then begin
      print, 'CREG not set, assuming center of image'
      creg = (sz[0]-1)/2
  endif

  if arg_present( SKY ) then sky = fltarr(sz[0], sz[1])
  if arg_present( SKYRMS ) then skyrms = fltarr(sz[0])

; Final image
  fimg = fltarr(sz[0], sz[1])

;
; Set region min/max

  if keyword_set( REG ) then begin
      szreg = size(reg, /dimensions)
      nreg = szreg[0]
;      regmin = min(reg[*,0]) > 0.
;      regmax = max(reg[*,1]) < float(sz[1]-1)
      regmin =  0.
      regmax = float(sz[1]-1)
      nwreg = reg    ; Create nwreg
  endif

;
; Set center of trace

  if keyword_set( TRACE ) then center = trace[long(creg)]

;;;;;;;;;;;;
; MAIN LOOP
  
  for q=0L,sz[0]-1 do begin 

      ; Set Region
      if keyword_set( REG ) then begin
          if keyword_set( TRACE ) then offset = trace[q] - center $
          else offset=0.
          for i=0,nreg-1 do begin
              for j=0,1 do begin   ; Restrain regions
                  nwreg[i,j] = regmax < (reg[i,j] + offset) > regmin 
              endfor
          endfor
      endif

      ; FIT
      if keyword_set( NOREJ ) then $ ;<NOREJ
        fit = x_fit(findgen(sz[1]), img[q,*], FITSTR=skyfstr, REG=nwreg) $
      else $  ; REJ
        fit = x_fitrej(findgen(sz[1]), img[q,*], FITSTR=skyfstr, REG=nwreg) 

      ; SUBTRACT THE SKY
      fimg[q,*] = img[q,*] - fit

      ; SAVE SKY
      if arg_present( SKY ) then sky[q,*] = fit
      if arg_present( SKYRMS ) then skyrms[q] = skyfstr.rms

  endfor

  delvarx, fit

  return, fimg
end
