;+ 
; NAME:
; wfccd_editobjstr
;    Version 1.0
;
; PURPOSE:
;    Finds position of science objects and seredinps in the image
;      and creates the object structure
;
; CALLING SEQUENCE:
;   
;   wfccd_editobjstr, wfccd, mask_id, [exp_id], /NOCHK
;
; INPUTS:
;   wfccd      - wfccd structure
;   mask_id    - Mask id value
;   exp_id     - Exposure id
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
;   wfccd_editobjstr, wfccd, mask_id, 0L
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro wfccd_editobjstr, wfccd, mask_id, exp_id, SILENT=silent, $
                    NOSLIT=noslit, NOOBJ=noobj

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_editobjstr, wfccd, mask_id, [exp_id], /SILENT, /NOSLIT'
    print, '         /NOOBJ [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( SKYOFF ) then skyoff = 2.

; Set exp
  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nexp)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp[0]

  if not keyword_set(SILENT) then print, 'wfccd_editobjstr: Reading...'

; Input the slit structure
  wfslit = xmrdfits(wfccd[exp].slit_fil,1,STRUCTYP='mslitstrct', /silent)

; Input the Obj structure
  objfil = 'Extract/Obj_'+strmid(wfccd[exp].img_root,3,3)+'.fits'
  wfobj = xmrdfits(objfil, 1, STRUCTYP='specobjstrct',/silent)

;  Read in Flux and Wave images
  flux = xmrdfits(wfccd[exp].img_final, /silent)
  wave = xmrdfits(wfccd[exp].img_final,2, /silent)

;  Check the Objects
  x_setobjgui, flux, wfobj, wfslit, wvimg=wave

; FITS
  ; WRITE
  if not keyword_set(NOOBJ) then mwrfits, wfobj, objfil, /create
  if not keyword_set(NOSLIT) then mwrfits, wfslit, wfccd[exp].slit_fil, /create
   
  ; ALL DONE
  if not keyword_set(SILENT) then print, 'wfccd_editobjstr: All Done!'
  return
end
