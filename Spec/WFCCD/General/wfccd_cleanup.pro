;+ 
; NAME:
; wfccd_cleanup   
;   Version 1.0
;
; PURPOSE:
;    Erases a number of files which are generally not needed
;      after reduction is complete
;
; CALLING SEQUENCE:
;   
;   wfccd_cleanup, wfccd, mask_id, flg
;
; INPUTS:
;   wfccd -- wfccd_strct defining the images of interest
;   mask_id -- Mask ID (e.g. 0L)
;   flg    -- Binary flag  (1=Flats, 2=Arcs, 4=Maps)
;
; RETURNS:
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
;   wfccd_cleanup, wfccd, 0L, 1
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   26-July-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_cleanup, wfccd, mask_id, flg

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'wfccd_cleanup, wfccd, mask_id, flg (v1.0)'
      return
  endif 
  
;  Optional Keywords
  
; FLATS
  if flg mod 2 EQ 1 then begin
      if mask_id LT 10 then flatfil = 'Flats/Flat_0'+$
        string(mask_id,format='(i1)')+'.fits' $
      else flatfil = 'Flats/Flat_'+string(mask_id,format='(i2)')+'.fits'
      a = findfile(flatfil+'*', count=nfil)
      if nfil NE 0 then spawn, '\rm '+flatfil+'*'
  endif

; ARCS
  if flg mod 4 GT 1 then begin
      arcs = where(wfccd.type EQ 'ARC' AND wfccd.flg_anly NE 0 AND $
                   wfccd.mask_id EQ mask_id, narc)
      for i=0L, narc-1 do begin
          ; I file
          arcifil = 'Arcs/ArcI_'+strmid(wfccd[arcs[i]].img_root,3)
          a = findfile(arcifil+'*', count=nfil)
          if nfil NE 0 then spawn, '\rm '+arcifil+'*'
          ; R file
          arcrfil = 'Arcs/ArcR_'+strmid(wfccd[arcs[i]].img_root,3)
          a = findfile(arcrfil+'*', count=nfil)
          if nfil NE 0 then begin
              spawn, '\rm '+arcrfil+'*'
              print, 'wfccd_cleanup:  Removing arcs'
          endif
      endfor
  endif

; MAPS
  if flg mod 8 GT 3 then begin
      stop
      if mask_id LT 10 then mapfil = 'Maps/Map_0'+$
        string(mask_id,format='(i1)')+'.fits' $
      else mapfil = 'Maps/Map_'+string(mask_id,format='(i2)')+'.fits'
      a = findfile(mapfil+'*', count=nfil)
      if nfil NE 0 then begin
          spawn, '\rm '+mapfil
          print, 'wfccd_cleanup:  Removing maps'
      endif
  endif

  return
end
