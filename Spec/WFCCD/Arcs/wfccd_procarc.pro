;+ 
; NAME:
; wfccd_procarc
;    Version 1.0
;
; PURPOSE:
;    Processes the Arc images
;
; CALLING SEQUENCE:
;   
;   wfccd_procarc, wfccd, WFAIMG=
;
; INPUTS:
;   wfccd    - WFCCD structure
;   mask_id     - Long defining the mask to process
;
; RETURNS:
;
; OUTPUTS:
;   wfaimg      -  WFCCD arc image (fits file)
;
; OPTIONAL KEYWORDS:
;  SLITSTR   - slit structure
;  WFASTR    - WFCCD arc structure
;  NOFITS   - Suppress fits output
;  OUTFIL   - Suppress fits output
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_procarc, wfccdstr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  MAIN DRIVER
pro wfccd_procarc, wfccd, mask_id, SVOV=svov, MAP=map, CLOBBER=clobber, $
                   FLAT=flat

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_procarc, wfccd, mask_id, /SVOV [v1.0]'
    return
  endif 

;  Optional Keywords
  a = findfile('Arcs/..', count=count)
  if count EQ 0 then file_mkdir, 'Arcs'

;  Find the Arc frames

  arcs = where(wfccd.type EQ 'ARC' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, narc)
  obj = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nobj)

; Check to see if Arcs exist
  if not keyword_set( CLOBBER ) then begin
      flgobj = intarr(nobj)
      for i=0L,nobj-1 do begin
          afil = wfccd[obj[i]].arc_fil
          if strlen(strtrim(afil,2)) LT 5 then flgobj[i] = 1 else begin
              a = findfile(afil, count=count)
              if count GT 0 then begin
                  print, 'wfccd_procarc: Arc file (', strtrim(afil,2), $
                    ') exists for obj ', $
                    strtrim(obj[i],2)
                  print, '                  -- Use /CLOBBER to reprocess'
              endif else flgobj[i]=1
          endelse
      endfor
      ; Keep obj with only those unprocessed files
      a = where( flgobj EQ 1, nobj)
      if nobj GT 0 then obj = obj[a] else return
  endif

;  OV subtract
  ovflt = where(wfccd[arcs].flg_ov EQ 0, nov)
  if nov NE 0 then wfccd_over, wfccd, arcs[ovflt], ORDR=4

; Flat
  if not keyword_set( FLAT ) then flat = xmrdfits(wfccd[arcs[0]].flat_fil,0,/silent)

;  Input MAP
  if not keyword_set( MAP ) then map = xmrdfits(wfccd[arcs[0]].map_fil, 0, /silent)

  if not keyword_set( SILENT ) then print, 'wfccd_procarc: Processing the Arcs'

;  Check for closest arc for each obj
  if narc GT 1 then begin
      flgarc = intarr(narc)
      flgobj = intarr(nobj)
      ; Find closest Arc to object frame
      for i=0L,nobj-1 do begin
          imn = wfccd_minut(wfccd, obj[i], arcs)
          flgarc[imn] = 1
          flgobj[i] = imn
      endfor
      ; Loop over the good arcs
      for i=0L,narc-1 do begin
          if flgarc[i] NE 0 then begin
              ; Read in OV image
              ovarc = xmrdfits(wfccd[arcs[i]].img_ov, /silent)
              ; Flatten
              ovarc = x_normspec(ovarc, flat)
              ; Undistort
              nwarc = x_rectify(ovarc, map, /silent)
              ; Output
              outfil = strtrim('Arcs/ArcR_'+strmid(wfccd[arcs[i]].img_root,3),2)
              mwrfits, nwarc, outfil, /create
              spawn, 'gzip -f '+outfil
              ; Update structure
              a = where(flgobj EQ i)
              wfccd[obj[a]].arc_fil = outfil
          endif
      endfor
  endif else begin
      ; Read in OV image
      ovarc = xmrdfits(wfccd[arcs[0]].img_ov, /silent)
      ; Flatten
      ovarc = x_normspec(ovarc, flat)
      ; Undistort
      nwarc = x_rectify(ovarc, map)
      ; Output
      outfil = 'Arcs/ArcR_'+strmid(wfccd[arcs[0]].img_root,3)
      mwrfits, nwarc, outfil, /create
      spawn, 'gzip -f '+outfil
      ; Update structure
      wfccd[obj].arc_fil = outfil
  endelse

  ; Delete OV files
  if not keyword_set( SVOV ) then wfccd_delov, wfccd, arcs
      
  return
end
