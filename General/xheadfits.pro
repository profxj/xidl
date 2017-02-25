;+ 
; NAME:
; xheadfits
;   Version 1.1
;
; PURPOSE:
;    Uses headfits to read the header of a fits file.  The only
;   advantage of xheadfits is that the file can be compressed and
;   yet the filename may be named without the gz extension.
;
; CALLING SEQUENCE:
;   head = xheadfits(fil, _EXTRA)
;
; INPUTS:
;   fil -- FITS Filename
;
; RETURNS:
;   head -- Header
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
;   head = xheadfits('spec.fits')
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function xheadfits, fil, _EXTRA=extra

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'head = xheadfits(fil, _EXTRA=) [V1.0]'
    return, -1
  endif 

;  Optional Keywords
;  if not keyword_set( extension ) then extension = 0L

; Filename
  if strlen(strtrim(fil,2)) EQ 0 then begin
      print, 'xmrdfits: File not set!'
      stop
      return, -1
  endif
  

  a = findfile(fil, count=na)
  if na EQ 0 then begin
      gz_fil = strtrim(fil,2)+'.gz'
      b = findfile(gz_fil, count=nb)
      if nb EQ 0 then begin
          print, 'xheadfits: Files ', fil, gz_fil, ' do not exist!'
          z_fil = strtrim(fil,2)+'.Z'
          c = findfile(z_fil, count=nc)
          if nc NE 0 then datfil = z_fil else return, -1 
      endif else datfil = gz_fil
  endif else datfil = fil

  ;; COMPRESS?
  len = strlen(fil)
;  if strmatch(strmid(fil,len-1),'z',/fold_case) $
;    then COMPRESS=1 else compress = 0

  return, headfits(datfil, COMPRESS=compress, _EXTRA=EXTRA)

end
