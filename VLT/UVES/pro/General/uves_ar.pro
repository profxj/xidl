;+ 
; NAME:
; uves_ar   
;     Version 1.1
;
; PURPOSE:
;   Reads in the first fits file in the directory with name
;   'uves_*.fits' and passes back the uves structure.
;
; CALLING SEQUENCE:
;   
;  uves = uves_ar(file)
;
; INPUTS:
;    [file] - Filename (default: first file in list ./uves_*fits*)
;
; RETURNS:
;    uves -  MIKE structure
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
;   uves = uves_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   uves_rslvall
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function uves_ar, file

; Optional Keywords
  if not keyword_set(file) then begin
      a = findfile('./uves_*fits*', count=count)
      if count EQ 0 then return, -1 
      print, 'Reading from '+a[0]
      file = a[0]
  endif

; Open and deal
  uves = xmrdfits(file,1,STRUCTYP='uvesstrct')
  ;; TRIM but allow for null
  uves.img_root = strtrim(uves.img_root,2)
  if max(strlen(uves.img_root)) EQ 0 then uves.img_root = ' '
  uves.obj_fil = strtrim(uves.obj_fil,2)
  if max(strlen(uves.obj_fil)) EQ 0 then uves.obj_fil = ' '
  uves.arc_fil = strtrim(uves.arc_fil,2)
  if max(strlen(uves.arc_fil)) EQ 0 then uves.arc_fil = ' '
;  uves.map_fil = strtrim(uves.map_fil,2)
;  if max(strlen(uves.map_fil)) EQ 0 then uves.map_fil = ' '
  uves.flat_fil = strtrim(uves.flat_fil,2)
  if max(strlen(uves.flat_fil)) EQ 0 then uves.flat_fil = ' '
  uves.Obj = strtrim(uves.Obj,2)
  if max(strlen(uves.Obj)) EQ 0 then uves.Obj = ' '
  uves.RA = strtrim(uves.RA,2)
  if max(strlen(uves.RA)) EQ 0 then uves.RA = ' '
  uves.ccd = strtrim(uves.ccd,2)
  if max(strlen(uves.ccd)) EQ 0 then uves.ccd = ' '
  uves.type = strtrim(uves.type,2)
  if max(strlen(uves.type)) EQ 0 then uves.type = ' '
  uves.img_final = strtrim(uves.img_final,2)
  if max(strlen(uves.img_final)) EQ 0 then uves.img_final = ' '
  uves.img_ov = strtrim(uves.img_ov,2)
  if max(strlen(uves.img_ov)) EQ 0 then uves.img_ov = ' '
;  uves.ccdspeed = strtrim(uves.ccdspeed,2)
;  if max(strlen(uves.ccdspeed)) EQ 0 then uves.ccdspeed = ' '

  ;; Resolve
  uves_rslvall

  return, uves

end
      
