;+
;
; NAME
;     copy_header.pro
;
; PURPOSE
;     Copies all of the keywords from a given header to a clean
;     header. The only keyword not copied is the 'EXTNAME'
;     keyword. This must be specified using the optional parameter
;     called "extname". 
;
;
; SYNTAX
;     newhdr = copy_header(hdr, [extname, comment=comment])
;
; INPUTS
;     hdr = the header which the user wishes to copy.
;     extname = a string giving the name of the extension name. The
;               header keyword 'EXTNAME' will be set equal to this
;               string.
;     comment = a comment string that will be written to the header to
;               comment the naming of the extension.
;
; KEYWORDS
;     None.
;
; OUTPUTS
;     A fresh clean header with all of the same keywords as the input
;     header (except that the 'DATE' keyword wil be
;     updated/added). This header will be ready to be written to a
;     file using mwrfits.pro.
;
; PROCEDURES CALLED 
;     fxhclean
;     fxhmake
;     sxaddpar
;     sxdelpar
;
; EXAMPLES
;
;
; COMMENTS
;
;
; HISTORY
;     Created February 4, 2003 by mcc.
;
;-

function copy_header, hdr, extname, comment=comment

; clean the input header stripping the header keywords which specify
; the data format.
  errmsg = ''
  fxhclean, hdr, err=errmsg

; restamp the date.
  get_date, date
  fxaddpar, hdr, 'DATE', date, 'Creation date'

; if the extname parameter is supplied, then add the extension name to
; the header. otherwise, remove any old extname from the header
; keywords.
  if n_params() gt 1 then begin
      if keyword_set(comment) then $
        sxaddpar, hdr, 'EXTNAME', extname, comment $
      else sxaddpar, hdr, 'EXTNAME', extname
  endif else sxdelpar, hdr, 'EXTNAME'

; output the modified header.
  return, hdr


end

