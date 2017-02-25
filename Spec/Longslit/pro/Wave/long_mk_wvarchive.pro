;+
; NAME:
;  long_mk_wvarchive
;  Version 1.1
;
; PURPOSE:  Generate an Archive file from a wave save file
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
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
;   Nov-2012  Written by JXP
;-
; long_mk_wvarchive, 'wave-r121019_0133.sav', 'lris_red_600_7500_red.sav'
;------------------------------------------------------------------------------
PRO LONG_MK_WVARCHIVE, new_arc, archive, ITAKE=itake, NOSTOP=nostop

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'long_addtoarchive, new_arc, archive, ITAKE= [v1.1]'
      return
  endif 

  
  calib_path = GETENV('LONGSLIT_DIR') + '/calib/linelists/'

  ;; Warning message
  print, 'long_addtoarchive:  Proceed only if: '
  print, 'long_addtoarchive:    1.  You know what you are doing'
  print, 'long_addtoarchive:    2.  You have notified JXP or JFH'
  print, 'long_addtoarchive:    3.  Your new arc (only takes one) is good'
  print, 'long_addtoarchive:    4.  Your new arc matches the npix of the archive'
  print, 'long_addtoarchive:  Continue as appropriate...'


  ;; Begin
  restore, new_arc
  if not keyword_set(NOSTOP) then stop

  sz = size(arc1d, /dimen)
  if n_elements(sz) GT 1 then begin
     if not keyword_set(ITAKE) then itake = lindgen(sz[1])
     arc1d = arc1d[*,itake]     ; Only takes the first element
     xfit = xfit[itake]
     stop
  endif

  ;; Open arcvhive
  archive_arc = arc1d
  calib = xfit

  ;; Writing
  save, archive_arc, calib, filen=calib_path+archive

  print, 'long_mk_wvarchive:  All done!'

  return
end
