;+ 
; NAME:
; hires_ar   
;     Version 1.1
;
; PURPOSE:
;   Reads in the first fits file in the directory with name
;   'hires_*.fits' and passes back the hires structure.
;
; CALLING SEQUENCE:
;   
;  hires = hires_ar(file)
;
; INPUTS:
;    [file] - Filename (default: first file in list ./hires_*fits*)
;
; RETURNS:
;    hires -  HIRES structure
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
;   hires = hires_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   hires_rslvall
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function hires_ar, file

; Optional Keywords
  if not keyword_set(file) then begin
      a = findfile('./hires_*fits*', count=count)
      if count EQ 0 then return, -1 
      print, 'Reading from '+a[0]
      file = a[0]
  endif

; Open and deal
  hires = xmrdfits(file,1,STRUCTYP='hiresstrct')
  ;; TRIM but allow for null
  hires.img_root = strtrim(hires.img_root,2)
  if max(strlen(hires.img_root)) EQ 0 then hires.img_root = ' '
  hires.obj_fil = strtrim(hires.obj_fil,2)
  if max(strlen(hires.obj_fil)) EQ 0 then hires.obj_fil = ' '
  hires.arc_fil = strtrim(hires.arc_fil,2)
  if max(strlen(hires.arc_fil)) EQ 0 then hires.arc_fil = ' '
;  hires.map_fil = strtrim(hires.map_fil,2)
;  if max(strlen(hires.map_fil)) EQ 0 then hires.map_fil = ' '
  hires.flat_fil = strtrim(hires.flat_fil,2)
  if max(strlen(hires.flat_fil)) EQ 0 then hires.flat_fil = ' '
  hires.Obj = strtrim(hires.Obj,2)
  if max(strlen(hires.Obj)) EQ 0 then hires.Obj = ' '
  hires.RA = strtrim(hires.RA,2)
  if max(strlen(hires.RA)) EQ 0 then hires.RA = ' '
  hires.ccd = strtrim(hires.ccd,2)
  if max(strlen(hires.ccd)) EQ 0 then hires.ccd = ' '
  hires.type = strtrim(hires.type,2)
  if max(strlen(hires.type)) EQ 0 then hires.type = ' '
  hires.img_final = strtrim(hires.img_final,2)
  if max(strlen(hires.img_final)) EQ 0 then hires.img_final = ' '
  hires.img_ov = strtrim(hires.img_ov,2)
  if max(strlen(hires.img_ov)) EQ 0 then hires.img_ov = ' '
;  hires.ccdspeed = strtrim(hires.ccdspeed,2)
;  if max(strlen(hires.ccdspeed)) EQ 0 then hires.ccdspeed = ' '

  ;; Resolve
  hires_rslvall

  return, hires

end
      
