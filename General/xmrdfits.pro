;+ 
; NAME:
; xmrdfits
;   Version 1.0
;
; PURPOSE:
;    Convert input to data whether it is a fits file or not
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
;  FSCALE      - Data is float
;
; OPTIONAL OUTPUTS:
;  HEAD        - Header
;
; COMMENTS:
;
; EXAMPLES:
;   dat = xmrdfits('spec.fits')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function xmrdfits, fil, extension, header, _EXTRA=extra

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'dat = xmrdfits(img, _EXTRA=) [V1.0]'
    return, -1
  endif 

;  Optional Keywords
  if not keyword_set( extension ) then extension = 0L

; Filename
  if strlen(fil) EQ 0 then begin
      print, 'xmrdfits: File not set!'
      stop
      return, -1
  endif
  

  a = findfile(fil, count=na)
  if na EQ 0 then begin
      gz_fil = fil+'.gz'
      b = findfile(gz_fil, count=nb)
      if nb EQ 0 then begin
          print, 'xmrdfits: Files ', fil, gz_fil, ' do not exist!'
          return, -1
      endif else datfil = gz_fil
  endif else datfil = fil

  return, mrdfits(datfil, extension, header, _EXTRA=EXTRA)

end
