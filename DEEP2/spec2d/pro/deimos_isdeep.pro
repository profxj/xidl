; file = optional parameter allowing user to pass name of file from
; which to grab the 'MaskDesign' table. this can either by a raw
; deimos frame or the bintab file. if not supplied, then the routine
; simply searches for the bintab file. obviously in creating the
; bintab file (within make)bintab_file.pro), one would need to pass a
; rawframe as this file parameter.

; isdeep = will be returned as 1 if the PROJNAME = 'DEEP2-1HS'

; maskname = if present as a second argument, then the maskname will
; be trimmed. the variable passed will be set to this trimmed string
; if trimming occurs or left alone if no trimming occurs.

pro deimos_isdeep, isdeep, maskname, file=file

; if a science frame or the bintab file isn't passed as the optional
; file parameter, then find the bintab file.
  if keyword_set(file) then $
    file = findfile(file, count=nfile) $
  else file = findfile('*bintab*.fits*', count=nfile)
  if nfile eq 0 then begin
      print, '(deimos_isdeep) ERROR: No bintab file found!!!'
      isdeep = 0
      return
  endif
; open the bintab file, get the extension names, and then close the
; file.
  file = file[0]
  fits_open, file, fcb
  extnames = strcompress(fcb.extname, /rem)
  fits_close, fcb
; find the 'MaskDesign' bin table.
  maskdex = where(extnames eq 'MaskDesign', maskcnt)
  if maskcnt eq 0 then begin
      print, '(deimos_isdeep) ERROR: No MaskDesign table found!!!'
      isdeep = 0
      return
  endif
  
  if getenv('IDLDEEPPHOTOMETRY_DIR') eq ""  then begin
    print, '(deimos_isdeep) ERROR (if DEEP data): No PCATs found!!!'
    isdeep = 0
    return
 endif

; extract the MaskDesign table.
  maskdex = maskdex[0]
  masktab = mrdfits(file, maskdex, /silent)
  project = strcompress(masktab.PROJNAME, /rem)
  if project eq 'DEEP2-1HS' then isdeep = 1 else isdeep = 0
; check if this is a KTRS mask. these masks will have a project name
; of DEEP2-1HS, but should be considere isdeep=0 masks.
  if n_params() gt 1 and size(maskname, /tname) eq 'STRING' then begin
      if strpos(strupcase(maskname), 'KTRS') ge 0 then isdeep = 0
  endif else begin
;      mask = (strsplit(file, '.', /extract))[0]
      hdr = headfits(file, /silent)
      mask = sxpar(hdr, 'SLMSKNAM')
      if strpos(strupcase(mask), 'KTRS') ge 0 then isdeep = 0
  endelse
; if a maskname has been passed, then trim it accordingly.
  if n_params() gt 1 and size(maskname, /tname) eq 'STRING' then begin
; check for the trailing .E's and .W's.
      maskname = strcompress(maskname, /rem)
      if stregex(maskname, '\.[EW]$') ne -1 then $
        if isdeep and (strlen(maskname) gt 4) then $
        maskname = strmid(maskname, 0, 4)
  endif

; return so that the isdeep variable returns set to either 0 or 1.
  return

end


