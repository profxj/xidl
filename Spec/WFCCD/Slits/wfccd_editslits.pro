;+ 
; NAME:
; wfccd_editslits
;    Version 1.0
;
; PURPOSE:
;    Parses the mask file for slit information
;
; CALLING SEQUENCE:
;   
;   wfccd_editslits, wfccd, mask_id, SLISTR=
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
;   wfccd_editslits, wfccd, mask_id, SLITSTR=
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  X_STSLTGUI
;
; REVISION HISTORY:
;   27-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro wfccd_editslits, wfccd, mask_id, flg, SLITSTR=slitstr, NOVERIFY=noverify, $
                    FLAT=flat, MAP=map, SVFLT=svflt, OUTFIL=outfil, SILENT=silent

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_editslits, wfccd, mask_id, [flg], SLITSTR=, NOVERIFY= [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( flg ) then flg = 0

; Obj and Flats
  flts = where(wfccd.type EQ 'FLT' AND wfccd.flg_anly NE 0 AND $
               wfccd.mask_id EQ mask_id, nflt)  
  obj = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nobj)

; Open the slit file
  slitstr = xmrdfits(wfccd[obj[0]].slit_fil, 1, structyp='slitstrct')

;  Read in the flat
  if not keyword_set( FLAT ) then $
    flat = xmrdfits(strtrim(wfccd[flts[0]].flat_fil,2), 0, /silent)

;  Read in the map
  map = xmrdfits(strtrim(wfccd[obj[0]].map_fil,2), 0, /silent)

; Undistort the flat 
  if not keyword_set( SILENT ) then $
    print, 'wfccd_editslits: Undistorting the flat...'
  nwflt = x_rectify(flat, map)

; Run the Gui
  x_stsltgui, nwflt, slitstr

;  Map the Slits into the Original Frame
  if not keyword_set( SILENT ) then $
    print, 'wfccd_editslits: Mapping the slits onto the original image...'
  imap = xmrdfits(strtrim(wfccd[obj[0]].map_fil,2), 1, /silent)
  x_origslit, slitstr, imap, /inverse

;  Output the Slits
  a = findfile('Slits/..', count=count)
  if count EQ 0 then file_mkdir, 'Slits'

  ; FITS
  outfil = wfccd[obj[0]].slit_fil
  mwrfits, slitstr, outfil, /create
   
  if not keyword_set(SILENT) then print, 'wfccd_setslits: All Done!'
  return
end
