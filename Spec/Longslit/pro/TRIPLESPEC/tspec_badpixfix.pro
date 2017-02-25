;+
; NAME:
;        TSPEC_BADPIXFIX
;
; PURPOSE:
;   Apply a bad pixel mask.  
;   Probably from the archived flat
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS
;
; SIDE EFFECTS:
;
;
; MODIFICATION HISTORY:
;  Written by JXP (August 2013) 
;       Version 1.0
;-
pro tspec_badpixfix, flux, invvar, BADPIXFIL=badpixfil

; Set the default xtalk factor to 170
  if not keyword_set(BADPIXFIL) then begin
     badpixfil = getenv('LONGSLIT_DIR')+'/pro/TRIPLESPEC/tspec-badpix-20120401_300s.fits'
     print, 'tspec_badpixfix: Using the default in the XIDL archive.'
  endif
  print, 'tspec_badpixfix: Using ', badpixfil
  
  ;; Read
  mask = xmrdfits(badpixfil,/sile)

  ;; 
  gdp = where(mask, complement=badp, ncomple=nbadp)
  sz = size(flux, /dimen)
  newfx = fltarr(sz[0], sz[1])
  newfx = flux[gdp]
  if keyword_set(invvar) and nbadp GT 0 then $
     invvar[badp] = -1.

  return
end

  
