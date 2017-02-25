;+ 
; NAME:
; imacsls_ar   
;     Version 1.1
;
; PURPOSE:
;   Reads in the first file in the directory with 'imacsls_*fits*'
;
; CALLING SEQUENCE:
;  imacsls = imacsls_ar(file)
;
; INPUTS:
;    [file] - Filename (default: first file in list ./imacsls_*fits*)
;
; RETURNS:
;    imacsls -  IMACS long slit structure
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
;   imacsls = imacsls_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Dec-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function imacsls_ar, file

; Optional Keywords
  if not keyword_set(file) then begin
      a = findfile('./imacsls_*fits*', count=count)
      if count EQ 0 then return, -1 
      print, 'Reading from '+a[0]
      file = a[0]
  endif

; Open and deal
  imacsls = xmrdfits(file,1,STRUCTYP='imacslsstrct')
  ;; TRIM but allow for null
  imacsls.img_root = strtrim(imacsls.img_root,2)
  if max(strlen(imacsls.img_root)) EQ 0 then imacsls.img_root = ' '
  imacsls.obj_fil = strtrim(imacsls.obj_fil,2)
  if max(strlen(imacsls.obj_fil)) EQ 0 then imacsls.obj_fil = ' '
  imacsls.arc_fil = strtrim(imacsls.arc_fil,2)
  if max(strlen(imacsls.arc_fil)) EQ 0 then imacsls.arc_fil = ' '
  imacsls.map_fil = strtrim(imacsls.map_fil,2)
  if max(strlen(imacsls.map_fil)) EQ 0 then imacsls.map_fil = ' '
  imacsls.flat_fil = strtrim(imacsls.flat_fil,2)
  if max(strlen(imacsls.flat_fil)) EQ 0 then imacsls.flat_fil = ' '
  imacsls.Obj = strtrim(imacsls.Obj,2)
  if max(strlen(imacsls.Obj)) EQ 0 then imacsls.Obj = ' '
  imacsls.RA = strtrim(imacsls.RA,2)
  if max(strlen(imacsls.RA)) EQ 0 then imacsls.RA = ' '
  imacsls.ccd = strtrim(imacsls.ccd,2)
  if max(strlen(imacsls.ccd)) EQ 0 then imacsls.ccd = ' '
  imacsls.type = strtrim(imacsls.type,2)
  if max(strlen(imacsls.type)) EQ 0 then imacsls.type = ' '
  imacsls.img_final = strtrim(imacsls.img_final,2)
  if max(strlen(imacsls.img_final)) EQ 0 then imacsls.img_final = ' '
  imacsls.img_ov = strtrim(imacsls.img_ov,2)
  if max(strlen(imacsls.img_ov)) EQ 0 then imacsls.img_ov = ' '
  imacsls.ccdspeed = strtrim(imacsls.ccdspeed,2)
  if max(strlen(imacsls.ccdspeed)) EQ 0 then imacsls.ccdspeed = ' '
  imacsls.grising = strtrim(imacsls.grising,2)
  if max(strlen(imacsls.grising)) EQ 0 then imacsls.grising = ' '

  ;; Resolve
;  imacsls_rslvall

  return, imacsls

end
      
