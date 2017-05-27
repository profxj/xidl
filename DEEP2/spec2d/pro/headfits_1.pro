function headfits_1, filename,  exten=exten
;+
; NAME:
;   headfits_1
; PURPOSE:
;   returns only first 36 lines of fits header.  This is useful for
;   DEIMOS headers which are very long, and where critical information
;   is in the beginning of the header
;
; CALLING SEQUENCE:
;   hdr = headfits_1(filename, exten=exten)
;
; INPUTS:
;   filename = (str) full filename
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;   exten  = (default 0) for primary header, or extension
;
; OUTPUTS:
;   hdr -- header from first block, 36 lines (2880 characters)
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;   md 27mar02
;
;----------------------------------------------------------------------

;position to correct header

if NOT keyword_set(exten) then exten_n = 0 else exten_n = exten
unit = fxposit( filename, exten_n, /READONLY, /SILENT)
 
if unit EQ -1 then begin 
      message = 'Unable to open file ' + filename 
      return,-1
endif
if eof(unit) then begin
      free_lun,unit
      message = 'Extension past EOF'
      return,-1
endif

;read first block
   hdr = bytarr(80, 36, /NOZERO)
   N_hdrblock = 1
   readu, unit, hdr
   header = string( hdr > 32b)

free_lun, unit
 
return, header

end
