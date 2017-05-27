function copy_keywords, header, keyword,  table_name
;+
; NAME:  
;   copy_keywords
;
; PURPOSE:
;   generates a new FITS binary table header and saves selected keywords
;
; CALLING SEQUENCE:
;   nhdr = copy_keywords(header, keyword,[table_name]) 
;
; INPUTS:
;   header  -- FITS header, original
;   keyword -- array of FITS keywords to be saved in new header
;
; OPTIONAL INPUTS:
;   table_name -- string name for binary table extension	
;
; KEYWORDS:
;
; OUTPUTS:
;   nhdr  -- new FITS header, with essential initial keywords, plus
;            copied keywords
;
; COMMENTS:
;   initializes new FITS header and adds DATE
;
; REVISION HISTORY:
;   MD 3jun2002
;----------------------------------------------------------------------

   extension = (n_params() gt 2) ? table_name : '' ;conditional statement
   fxbhmake, nhdr, 1,  extension, /initialize, /date

   for i= 0, n_elements(keyword)-1 do begin
      value = fxpar(header, keyword[i])
      sxaddpar, nhdr, keyword[i], value
   endfor

return, nhdr
end












