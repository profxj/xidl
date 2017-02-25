;+ 
; NAME:
; xmrdfits
;   Version 1.1
;
; PURPOSE:
;    Uses mrdfits to read the header of a fits file.  The only
;   advantage of xheadfits is that the file can be compressed and
;   yet the filename may be named without the gz extension.
;
; CALLING SEQUENCE:
;   
;   dat = xmrdfits(fil, [extension, header], _EXTRA)
;
; INPUTS:
;   fil -- FITS Filename
;   [extension] -- Extension of the data array
;
; RETURNS:
;   dat  - Data array
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  header -- Header of the FITS extension
;
; COMMENTS:
;
; EXAMPLES:
;   dat = xmrdfits('spec.fits')
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Sep-2002 Written by JXP
;   08-Jun-2011 Use IDL's file_search() instead of findfile()
;               to be MUCH faster, KLC
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function xmrdfits, fil, extension, header, _EXTRA=extra, FLG_FIL=FLG_FIL

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'dat = xmrdfits(fil, extension, header, _EXTRA=, FLG_FIL=) [V1.1]'
    return, -1
  endif 

;  Optional Keywords
  if not keyword_set( extension ) then extension = 0L

  FLG_FIL = 1L

; Filename
  if strlen(fil) EQ 0 then begin
      print, 'xmrdfits: File not set!'
      FLG_FIL = 0L
      stop
      return, -1
  endif
  

  a = file_search(fil, count=na)
  if na EQ 0 then begin
      gz_fil = strtrim(fil,2)+'.gz'
      b = file_search(gz_fil, count=nb)
      if nb EQ 0 then begin
          print, 'xmrdfits: Files ', fil, gz_fil, ' do not exist!'
          FLG_FIL = 0L
          ;; Quicker to do 2nd search only as needed
          z_fil = strtrim(fil,2)+'.Z'
          c = file_search(z_fil, count=nc)
          if nc NE 0 then datfil = z_fil else return, -1
      endif else datfil = gz_fil
  endif else datfil = fil

  return, mrdfits(datfil, extension, header, _EXTRA=EXTRA)

end
