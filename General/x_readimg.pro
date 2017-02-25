;+ 
; NAME:
; x_readimg   
;   Version 1.1
;
; PURPOSE:
;    Convert input to data whether it is a fits file or an image array
;
; CALLING SEQUENCE:
;   
;   dat = x_readimg(img)
;
; INPUTS:
;   img       - Fits file or data
;
; RETURNS:
;   dat       - Data in fits file or the input data
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /FSCALE      - Data is float
;  /DSCALE      - Data is double
;
; OPTIONAL OUTPUTS:
;  HEAD=       - Header
;
; COMMENTS:
;
; EXAMPLES:
;   dat = x_readimg('spec.fits')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_readimg, img, FSCALE=fscale, HEAD=head, DSCALE=dscale

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'dat = x_readimg(img, /fscale, /dscale, HEAD=) [V1.1]'
    return, -1
  endif 

;  Optional Keywords

  if keyword_set(DSCALE) and keyword_set(FSCALE) then fscale=0

; Test if it is a file or image

  if size(img, /type) EQ 7 then begin
      a = findfile(img+'*', count=count)
      if count EQ 0 then begin
          print, img+' does not exist'
          return, -1
      endif
      ; Read in ydat 
      dat = xmrdfits(img, 0, head, fscale=fscale, dscale=dscale, /silent)
  endif else begin
      if keyword_set( FSCALE ) OR keyword_set( DSCALE ) then begin
          if keyword_set( FSCALE ) then dat = float(img)
          if keyword_set( DSCALE ) then dat = double(img)
      endif else dat = img
  endelse
      
  return, dat
end
