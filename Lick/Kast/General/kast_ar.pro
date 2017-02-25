;+ 
; NAME:
; kast_ar   
;     Version 1.1
;
; PURPOSE:
;   Reads the Kast IDL structure from fits file into memory.  The
;  program grabs the first file in the directory with kast_*fits
;
; CALLING SEQUENCE:
;   
;  kast = kast_ar(file)
;
; INPUTS:
;    [file] - Filename (default: first file in list ./kast_*fits*)
;
; RETURNS:
;    kast -  Kast IDL structure
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
;   kast = kast_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function kast_ar, file

; Optional Keywords
  if not keyword_set(file) then begin
      a = findfile('./kast_*fits*', count=count)
      if count EQ 0 then return, -1 
      print, 'Reading from '+a[0]
      file = a[0]
  endif

; Open and deal
  kast = xmrdfits(file,1,STRUCTYP='kaststrct')
  ;; TRIM but allow for null
  kast.grising = strtrim(kast.grising,2)
  if max(strlen(kast.grising)) EQ 0 then kast.grising = ' '
  kast.obj_fil = strtrim(kast.obj_fil,2)
  if max(strlen(kast.img_root)) EQ 0 then kast.img_root = ' '
  kast.obj_fil = strtrim(kast.obj_fil,2)
  if max(strlen(kast.obj_fil)) EQ 0 then kast.obj_fil = ' '
  kast.arc_fil = strtrim(kast.arc_fil,2)
  if max(strlen(kast.arc_fil)) EQ 0 then kast.arc_fil = ' '
  kast.map_fil = strtrim(kast.map_fil,2)
  if max(strlen(kast.map_fil)) EQ 0 then kast.map_fil = ' '
  kast.flat_fil = strtrim(kast.flat_fil,2)
  if max(strlen(kast.flat_fil)) EQ 0 then kast.flat_fil = ' '
  kast.Obj = strtrim(kast.Obj,2)
  if max(strlen(kast.Obj)) EQ 0 then kast.Obj = ' '
  kast.TEL = strtrim(kast.TEL,2)
  if max(strlen(kast.TEL)) EQ 0 then kast.TEL = ' '
  kast.RA = strtrim(kast.RA,2)
  if max(strlen(kast.RA)) EQ 0 then kast.RA = ' '
  kast.ccd = strtrim(kast.ccd,2)
  if max(strlen(kast.ccd)) EQ 0 then kast.ccd = ' '
  kast.type = strtrim(kast.type,2)
  if max(strlen(kast.type)) EQ 0 then kast.type = ' '
  kast.img_final = strtrim(kast.img_final,2)
  if max(strlen(kast.img_final)) EQ 0 then kast.img_final = ' '
  kast.img_ov = strtrim(kast.img_ov,2)
  if max(strlen(kast.img_ov)) EQ 0 then kast.img_ov = ' '
  kast.ccdspeed = strtrim(kast.ccdspeed,2)
  if max(strlen(kast.ccdspeed)) EQ 0 then kast.ccdspeed = ' '
  return, kast

end
      
