;+ 
; NAME:
; esi_ar   
;     Version 1.1
;
; PURPOSE:
;   Reads in the first file in the directory with 'esi*fits'
;
; CALLING SEQUENCE:
;   
;  esi = esi_ar(file)
;
; INPUTS:
;    [file] - Filename (default: first file in list ./esi_*fits*)
;
; RETURNS:
;    esi -  ESI structure
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
;   esi = esi_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function esi_ar, file

; Optional Keywords
  if not keyword_set(file) then begin
      a = findfile('./esi_*fits*', count=count)
      if count EQ 0 then return, -1 
      print, 'Reading from '+a[0]
      file = a[0]
  endif

; Open and deal
  esi = xmrdfits(file,1,STRUCTYP='esistrct')
  ;; TRIM but allow for null
  esi.img_root = strtrim(esi.img_root,2)
  if max(strlen(esi.img_root)) EQ 0 then esi.img_root = ' '
  esi.obj_fil = strtrim(esi.obj_fil,2)
  if max(strlen(esi.obj_fil)) EQ 0 then esi.obj_fil = ' '
  esi.arc_fil = strtrim(esi.arc_fil,2)
  if max(strlen(esi.arc_fil)) EQ 0 then esi.arc_fil = ' '
  esi.map_fil = strtrim(esi.map_fil,2)
  if max(strlen(esi.map_fil)) EQ 0 then esi.map_fil = ' '
  esi.flat_fil = strtrim(esi.flat_fil,2)
  if max(strlen(esi.flat_fil)) EQ 0 then esi.flat_fil = ' '
  esi.Obj = strtrim(esi.Obj,2)
  if max(strlen(esi.Obj)) EQ 0 then esi.Obj = ' '
  esi.TEL = strtrim(esi.TEL,2)
  if max(strlen(esi.TEL)) EQ 0 then esi.TEL = ' '
  esi.RA = strtrim(esi.RA,2)
  if max(strlen(esi.RA)) EQ 0 then esi.RA = ' '
  esi.ccd = strtrim(esi.ccd,2)
  if max(strlen(esi.ccd)) EQ 0 then esi.ccd = ' '
  esi.type = strtrim(esi.type,2)
  if max(strlen(esi.type)) EQ 0 then esi.type = ' '
  esi.img_final = strtrim(esi.img_final,2)
  if max(strlen(esi.img_final)) EQ 0 then esi.img_final = ' '
  esi.img_ov = strtrim(esi.img_ov,2)
  if max(strlen(esi.img_ov)) EQ 0 then esi.img_ov = ' '
  esi.ccdspeed = strtrim(esi.ccdspeed,2)
  if max(strlen(esi.ccdspeed)) EQ 0 then esi.ccdspeed = ' '
  return, esi

end
      
