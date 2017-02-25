;+ 
; NAME:
; wfccd_setslits
;    Version 1.0
;
; PURPOSE:
;    Parses the mask file for slit information
;
; CALLING SEQUENCE:
;   
;   wfccd_setslits, wfccd, mask_id, SLISTR=
;
; INPUTS:
;   wfccd      - wfccd structure
;   mask_id    - Mask id value
;
; RETURNS:
;
; OUTPUTS:
;   slitstr      -  Creates a structure of slits
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_setslits, wfccd, mask_id, SLITSTR=
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  X_SETSLITS
;  X_STSLTGUI
;  X_RECTIFY
;  X_ORIGSLIT
;
; REVISION HISTORY:
;   19-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro wfccd_setslits, wfccd, mask_id, flg, SLITSTR=slitstr, NOVERIFY=noverify, $
                    FLAT=flat, MAP=map, SVFLT=svflt, OUTFIL=outfil, $
                    SILENT=silent, DEBUG=debug, THETA=theta, SHIFT=shift, $
                    NOFIT=nofit

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_setslits, wfccd, mask_id, [flg], SLITSTR=, NOVERIFY=, /DEBUG [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( flg ) then flg = 0

;  Find the files
  flts = where(wfccd.type EQ 'FLT' AND wfccd.flg_anly NE 0 AND $
               wfccd.mask_id EQ mask_id, nflt)  
  obj = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nobj)

; Parse the slit file
  slitstr = wfccd_prsslits(wfccd[obj[0]].msk_fil, flg)
  slitstr.field = strtrim(wfccd[obj[0]].Obj,2)

;  Read in the flat
  if not keyword_set( FLAT ) then $
    flat = xmrdfits(strtrim(wfccd[flts[0]].flat_fil,2), 0, /silent)

;  Read in the map
  if not keyword_set( MAP ) then $
    map = xmrdfits(strtrim(wfccd[obj[0]].map_fil,2), 0, /silent)

; Undistort the flat 
  if not keyword_set( SILENT ) then $
    print, 'wfccd_setslits: Undistorting the flat...'
  nwflt = x_rectify(flat, map)

;  Slits
  if not keyword_set( SILENT ) then $
    print, 'wfccd_setslits: Setting the slits...'
  if(wfccd[flts[0]].rotated) then begin
      splog, 'rotated slit, will not match flat!!'
      x_setslits, nwflt, slitstr, YSTRT=wfccd[flts[0]].ystrt, $
        DEBUG=debug, /nofit, theta=wfccd[flts[0]].theta, $
        shift=wfccd[flts[0]].shift
      if not keyword_set( NOVERIFY ) then begin
          splog, 'comparing to the exposure'
          img=xmrdfits(strtrim(wfccd[obj[0]].rootpth+ $
                               wfccd[obj[0]].img_root,2), 0, /silent)
          img=img[0:2047,0:2046]
          nx=(size(img,/dim))[0]
          img = x_rectify(img, map)
          x_stsltgui, img, slitstr
      endif
  endif else begin
      x_setslits, nwflt, slitstr, YSTRT=wfccd[flts[0]].ystrt, $
        DEBUG=debug
;  Check the slits
      if not keyword_set( NOVERIFY ) then x_stsltgui, nwflt, slitstr
  endelse


;  Map the Slits into the Original Frame
  if not keyword_set( SILENT ) then $
    print, 'wfccd_setslits: Mapping the slits onto the original image...'
  imap = xmrdfits(strtrim(wfccd[obj[0]].map_fil,2), 1, /silent)
  x_origslit, slitstr, imap, /inverse

;  Output the Slits
  a = findfile('Slits/..', count=count)
  if count EQ 0 then file_mkdir, 'Slits'

  if not keyword_set( OUTFIL ) then begin
      if mask_id LT 10 then outfil = 'Slits/Slits_0'+$
        string(mask_id,format='(i1)')+'.fits' $
      else outfil = 'Slits/Slits_'+string(mask_id,format='(i2)')+'.fits'
  endif
  ; Structure
  allimg = where(wfccd.mask_id EQ mask_id AND wfccd.flg_anly NE 0)
  wfccd[allimg].slit_fil = outfil
  ; FITS
  mwrfits, slitstr, outfil, /create
  spawn, 'gzip -f '+outfil
   
  if not keyword_set(SILENT) then print, 'wfccd_setslits: All Done!'
  return
end
