;+ 
; NAME:
; hires_makee_archive
;     Version 1.0
;
; PURPOSE:
;    Generate XIDL/HIRedux template files from the MAKEE database
;
; CALLING SEQUENCE:
;   
; INPUTS:
;              structure to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Oct-2012 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_makee_archive, MAKEEpath

;expects col=0 for oldest, col=1 for uv, col=2 for red (all pre
;mosaic)

  print, 'hires_makee_archive: Continue only if you know what you are doing!'
  stop

  fits_files = file_search(MAKEEpath+'ad*.fits', count=nfiles)

  ;; Loop
  for qq=0L,nfiles-1 do begin
     ;; Outfil
     ipos1 = strpos(fits_files[qq],'/', /reverse_search)
     ipos2 = strpos(fits_files[qq],'.fits')
     outfil = strmid(fits_files[qq], ipos1, ipos2-ipos1-1)+'_template.idl'
     outfil = getenv('HIRES_CALIBS')+'/ARCS/MAKEE/'+outfil

     hires_mk_makee_template, fits_files[qq], outfil
  endfor

  return

end
