;+ 
; NAME:
; wfccd_ar   
;     Version 1.0
;
; PURPOSE:
;   Reads in the first file in the directory with 'wfccd*fits'
;
; CALLING SEQUENCE:
;   
;  wfccd = wfccd_ar()
;
; INPUTS:
;
; RETURNS:
;    wfccd -  WFCCD structure
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
;   wfccd = wfccd_ar()
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function wfccd_ar

;  Optional Keywords
  
; Find the file
  a = findfile('./wfccd*fits', count=count)
  if count EQ 0 then return, -1 $
  else begin
      print, 'Reading from '+a[0]
      wfccd = mrdfits(a[0],1,STRUCTYP='wfccdstrct')
      ;; TRIM but allow for null
      wfccd.img_root = strtrim(wfccd.img_root,2)
      if max(strlen(wfccd.img_root)) EQ 0 then wfccd.img_root = ' '
      wfccd.obj_fil = strtrim(wfccd.obj_fil,2)
      if max(strlen(wfccd.obj_fil)) EQ 0 then wfccd.obj_fil = ' '
      wfccd.arc_fil = strtrim(wfccd.arc_fil,2)
      if max(strlen(wfccd.arc_fil)) EQ 0 then wfccd.arc_fil = ' '
      wfccd.slit_fil = strtrim(wfccd.slit_fil,2)
      if max(strlen(wfccd.slit_fil)) EQ 0 then wfccd.slit_fil = ' '
      wfccd.msk_fil = strtrim(wfccd.msk_fil,2)
      if max(strlen(wfccd.msk_fil)) EQ 0 then wfccd.msk_fil = ' '
      wfccd.map_fil = strtrim(wfccd.map_fil,2)
      if max(strlen(wfccd.map_fil)) EQ 0 then wfccd.map_fil = ' '
      wfccd.flat_fil = strtrim(wfccd.flat_fil,2)
      if max(strlen(wfccd.flat_fil)) EQ 0 then wfccd.flat_fil = ' '
      wfccd.filter = strtrim(wfccd.filter,2)
      if max(strlen(wfccd.filter)) EQ 0 then wfccd.filter = ' '
      wfccd.Obj = strtrim(wfccd.Obj,2)
      if max(strlen(wfccd.Obj)) EQ 0 then wfccd.Obj = ' '
      wfccd.TEL = strtrim(wfccd.TEL,2)
      if max(strlen(wfccd.TEL)) EQ 0 then wfccd.TEL = ' '
      wfccd.img_final = strtrim(wfccd.img_final,2)
      if max(strlen(wfccd.img_final)) EQ 0 then wfccd.img_final = ' '
      wfccd.img_ov = strtrim(wfccd.img_ov,2)
      if max(strlen(wfccd.img_ov)) EQ 0 then wfccd.img_ov = ' '
      wfccd.img_msk= strtrim(wfccd.img_msk,2)
      if max(strlen(wfccd.img_msk)) EQ 0 then wfccd.img_msk = ' '
      return, wfccd
  endelse

end
      
