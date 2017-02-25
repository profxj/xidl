;+ 
; NAME:
; mike_ar   
;     Version 1.1
;
; PURPOSE:
;   Reads in the first fits file in the directory with name
;   'mike_*.fits' and passes back the mike structure.
;
; CALLING SEQUENCE:
;   
;  mike = mike_ar(file)
;
; INPUTS:
;    [file] - Filename (default: first file in list ./mike_*fits*)
;
; RETURNS:
;    mike -  MIKE structure
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
;   mike = mike_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   mike_rslvall
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mike_ar, file

; Optional Keywords
  if not keyword_set(file) then begin
      a = findfile('./mike_*fits*', count=count)
      if count EQ 0 then return, -1 
      print, 'Reading from '+a[0]
      file = a[0]
  endif

; Open and deal
  mike = xmrdfits(file,1,STRUCTYP='mikestrct')
  ;; TRIM but allow for null
  mike.rootpth = strtrim(mike.rootpth,2)
  if max(strlen(mike.rootpth)) EQ 0 then mike.rootpth = ' '
  mike.img_root = strtrim(mike.img_root,2)
  if max(strlen(mike.img_root)) EQ 0 then mike.img_root = ' '
  mike.obj_fil = strtrim(mike.obj_fil,2)
  if max(strlen(mike.obj_fil)) EQ 0 then mike.obj_fil = ' '
  mike.arc_fil = strtrim(mike.arc_fil,2)
  if max(strlen(mike.arc_fil)) EQ 0 then mike.arc_fil = ' '
  mike.map_fil = strtrim(mike.map_fil,2)
  if max(strlen(mike.map_fil)) EQ 0 then mike.map_fil = ' '
  mike.flat_fil = strtrim(mike.flat_fil,2)
  if max(strlen(mike.flat_fil)) EQ 0 then mike.flat_fil = ' '
  mike.Obj = strtrim(mike.Obj,2)
  if max(strlen(mike.Obj)) EQ 0 then mike.Obj = ' '
  mike.RA = strtrim(mike.RA,2)
  if max(strlen(mike.RA)) EQ 0 then mike.RA = ' '
  mike.ccd = strtrim(mike.ccd,2)
  if max(strlen(mike.ccd)) EQ 0 then mike.ccd = ' '
  mike.type = strtrim(mike.type,2)
  if max(strlen(mike.type)) EQ 0 then mike.type = ' '
  mike.img_final = strtrim(mike.img_final,2)
  if max(strlen(mike.img_final)) EQ 0 then mike.img_final = ' '
  mike.img_ov = strtrim(mike.img_ov,2)
  if max(strlen(mike.img_ov)) EQ 0 then mike.img_ov = ' '
  mike.ccdspeed = strtrim(mike.ccdspeed,2)
  if max(strlen(mike.ccdspeed)) EQ 0 then mike.ccdspeed = ' '

  ;; Resolve
  mike_rslvall

  return, mike

end
      
