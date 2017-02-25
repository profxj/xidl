;+
; NAME:
;  long_addtoarchive
;  Version 1.1
;
; PURPOSE:
;
; CALLING SEQUENCE:
;  LONG_REIDENTIFY, arc_obj, lines, wstruct
;
; INPUTS:
;
; OPTIONAL INPUTS:
; ITAKE=  -- Index of the arc to add.  Default=0L
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
;
; REVISION HISTORY:
;   June-2009  Written by JXP
;-
; long_addtoarchive, 'wave-b1029.sav', 'kast_600_4310.sav'
;------------------------------------------------------------------------------
PRO LONG_ADDTOARCHIVE, new_arc, archive, ITAKE=itake, NOSTOP=nostop

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'long_addtoarchive, new_arc, archive, ITAKE= [v1.1]'
      return
  endif 

  if not keyword_set(ITAKE) then itake = 0L
  
  calib_path = GETENV('LONGSLIT_DIR') + '/calib/linelists/'

  ;; Warning message
  print, 'long_addtoarchive:  Proceed only if: '
  print, 'long_addtoarchive:    1.  You know what you are doing'
  print, 'long_addtoarchive:    2.  You have notified JXP or JFH'
  print, 'long_addtoarchive:    3.  Your new arc (only takes one) is good'
  print, 'long_addtoarchive:    4.  Your new arc matches the npix of the archive'
  print, 'long_addtoarchive:  Continue as appropriate...'
  if not keyword_set(NOSTOP) then stop


  ;; Begin
  restore, new_arc

  sz = size(arc1d, /dimen)
  if n_elements(sz) GT 0 then begin
     arc1d = arc1d[*,itake]         ; Only takes the first element
     xfit = xfit[itake]
  endif

  ;; Open arcvhive
  restore, calib_path+archive

  
  sza = size(archive_arc, /dimen)
  if sza[0] NE sz[0] then begin
      print, 'long_addtoarchive:  Your arc size does not match the archive!'
      return
  endif
  if n_elements(sza) GT 1 then nold = sza[1] else nold = 1

  archive_arc = reform([archive_arc[*], arc1d], [sza[0], nold+1])
  tmp = calib[0]
  copy_struct, xfit, tmp
  calib = [calib, tmp]

  ;; Writing
  save, archive_arc, calib, filen=calib_path+archive

  print, 'long_addtoarchive:  All done!'

  return
end
