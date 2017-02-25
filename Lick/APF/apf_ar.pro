;+ 
; NAME:
; apf_ar   
;     Version 1.1
;
; PURPOSE:
;   Reads in the first fits file in the directory with name
;   'apf_*.fits' and passes back the apf structure.
;
; CALLING SEQUENCE:
;   
;  apf = apf_ar(file)
;
; INPUTS:
;    [file] - Filename (default: first file in list ./apf_*fits*)
;
; RETURNS:
;    apf -  HIRES structure
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
;   apf = apf_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   apf_rslvall
;
; REVISION HISTORY:
;   18-Sep-2014 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function apf_ar, file

; Optional Keywords
  if not keyword_set(file) then begin
      a = findfile('./apf_*fits*', count=count)
      if count EQ 0 then return, -1 
      print, 'Reading from '+a[0]
      file = a[0]
  endif

; Open and deal
  apf = xmrdfits(file,1,STRUCTYP='apfstrct')
  ;; TRIM but allow for null
  apf.img_root = strtrim(apf.img_root,2)
  if max(strlen(apf.img_root)) EQ 0 then apf.img_root = ' '
  apfbj_fil = strtrim(apf.obj_fil,2)
  if max(strlen(apf.obj_fil)) EQ 0 then apf.obj_fil = ' '
  apf.arc_fil = strtrim(apf.arc_fil,2)
  if max(strlen(apf.arc_fil)) EQ 0 then apf.arc_fil = ' '
;  apf.map_fil = strtrim(apf.map_fil,2)
;  if max(strlen(apf.map_fil)) EQ 0 then apf.map_fil = ' '
  apf.flat_fil = strtrim(apf.flat_fil,2)
  if max(strlen(apf.flat_fil)) EQ 0 then apf.flat_fil = ' '
  apf.Obj = strtrim(apf.Obj,2)
  if max(strlen(apf.Obj)) EQ 0 then apf.Obj = ' '
  apf.RA = strtrim(apf.RA,2)
  if max(strlen(apf.RA)) EQ 0 then apf.RA = ' '
  apf.ccd = strtrim(apf.ccd,2)
  if max(strlen(apf.ccd)) EQ 0 then apf.ccd = ' '
  apf.type = strtrim(apf.type,2)
  if max(strlen(apf.type)) EQ 0 then apf.type = ' '
  apf.img_final = strtrim(apf.img_final,2)
  if max(strlen(apf.img_final)) EQ 0 then apf.img_final = ' '
  apf.img_ov = strtrim(apf.img_ov,2)
  if max(strlen(apf.img_ov)) EQ 0 then apf.img_ov = ' '
;  apf.ccdspeed = strtrim(apf.ccdspeed,2)
;  if max(strlen(apf.ccdspeed)) EQ 0 then apf.ccdspeed = ' '

  ;; Resolve
  apf_rslvall

  return, apf

end
      
