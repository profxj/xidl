;+ 
; NAME:
; hamspec_ar   
;     Version 1.1
;
; PURPOSE:
;   Reads in the first fits file in the directory with name
;   'hamspec_*.fits' and passes back the hamspec structure.
;
; CALLING SEQUENCE:
;   
;  hamspec = hamspec_ar(file)
;
; INPUTS:
;    [file] - Filename (default: first file in list ./hamspec_*fits*)
;
; RETURNS:
;    hamspec -  HIRES structure
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
;   hamspec = hamspec_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   hamspec_rslvall
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;   03-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function hamspec_ar, file

; Optional Keywords
  if not keyword_set(file) then begin
      a = findfile('./hamspec_*fits*', count=count)
      if count EQ 0 then return, -1 
      print, 'Reading from '+a[0]
      file = a[0]
  endif

; Open and deal
  hamspec = xmrdfits(file,1,STRUCTYP='hamspecstrct')
  ;; TRIM but allow for null
  hamspec.img_root = strtrim(hamspec.img_root,2)
  if max(strlen(hamspec.img_root)) EQ 0 then hamspec.img_root = ' '
  hamspec.obj_fil = strtrim(hamspec.obj_fil,2)
  if max(strlen(hamspec.obj_fil)) EQ 0 then hamspec.obj_fil = ' '
  hamspec.arc_fil = strtrim(hamspec.arc_fil,2)
  if max(strlen(hamspec.arc_fil)) EQ 0 then hamspec.arc_fil = ' '
;  hamspec.map_fil = strtrim(hamspec.map_fil,2)
;  if max(strlen(hamspec.map_fil)) EQ 0 then hamspec.map_fil = ' '
  hamspec.flat_fil = strtrim(hamspec.flat_fil,2)
  if max(strlen(hamspec.flat_fil)) EQ 0 then hamspec.flat_fil = ' '
  hamspec.pflat_fil = strtrim(hamspec.pflat_fil,2)
  if max(strlen(hamspec.pflat_fil)) EQ 0 then hamspec.pflat_fil = ' '
  hamspec.Obj = strtrim(hamspec.Obj,2)
  if max(strlen(hamspec.Obj)) EQ 0 then hamspec.Obj = ' '
  hamspec.RA = strtrim(hamspec.RA,2)
  if max(strlen(hamspec.RA)) EQ 0 then hamspec.RA = ' '
  hamspec.ccd = strtrim(hamspec.ccd,2)
  if max(strlen(hamspec.ccd)) EQ 0 then hamspec.ccd = ' '
  hamspec.type = strtrim(hamspec.type,2)
  if max(strlen(hamspec.type)) EQ 0 then hamspec.type = ' '
  hamspec.img_final = strtrim(hamspec.img_final,2)
  if max(strlen(hamspec.img_final)) EQ 0 then hamspec.img_final = ' '
  hamspec.img_ov = strtrim(hamspec.img_ov,2)
  if max(strlen(hamspec.img_ov)) EQ 0 then hamspec.img_ov = ' '
;  hamspec.ccdspeed = strtrim(hamspec.ccdspeed,2)
;  if max(strlen(hamspec.ccdspeed)) EQ 0 then hamspec.ccdspeed = ' '

  ;; Resolve
  hamspec_rslvall

  return, hamspec

end
      
