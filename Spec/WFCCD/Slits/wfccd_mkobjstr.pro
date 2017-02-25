;+ 
; NAME:
; wfccd_mkobjstr
;    Version 1.0
;
; PURPOSE:
;    Finds position of science objects and seredinps in the image
;      and creates the object structure
;
; CALLING SEQUENCE:
;   
;   wfccd_mkobjstr, wfccd, mask_id, [exp_id], /NOCHK
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
;   wfccd_mkobjstr, wfccd, mask_id, 0L
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Apr-2002 Written by JXP
;   12-May-2002 Revisions for sky edges
;-
;------------------------------------------------------------------------------

pro wfccd_mkobjstr, wfccd, mask_id, exp_id, NOCHK=nochk, SILENT=silent, $
                    NOSLIT=noslit, NOOBJ=noobj, NOSKYSET=noskyset, $ 
                    SKYOFF=skyoff

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_mkobjstr, wfccd, mask_id, [exp_id], /NOCHK, /SILENT, /NOSLIT'
    print, '         /NOOBJ, /NOSKYSET, SKYOFF= [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( SKYOFF ) then skyoff = 2.

  wfccd.img_final=strtrim(wfccd.img_final,2)

; Set exp
  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nexp)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp[0]

; Input the slit structure
  wfslit = xmrdfits(wfccd[exp].slit_fil,1,STRUCTYP='mslitstrct',/silent)

;  Read in Flux and Wave images
  splog, 'reading from '+wfccd[exp].img_final
  flux = xmrdfits(wfccd[exp].img_final, /silent)
  wave = xmrdfits(wfccd[exp].img_final,2, /silent)

; Check for objects Automatically and Create Obj Structure
  x_fndslitobj, flux, wave, wfslit, wfobj, WVGUESS=5720.
  wfobj.spec2d_fil = wfccd[exp].img_final
  wfobj.slit_fil = wfccd[exp].slit_fil
  wfobj.exp = wfccd[exp].exp

; Set yedg_sky
  if not keyword_set( NOSKYSET ) then begin
      if wfslit[0].yedg_sky[0,0] NE 0. and keyword_set(NOCHK) eq 0 then begin
          print, 'wfccd_mkobjstr: Are you sure you want to change the sky? (1/0)'
          ans = x_guinum(2, TITLE='Change sky? (1/0)')
      endif else ans = 1L
      if ans EQ 1L then begin
          wfslit.yedg_sky[*,0] = wfslit.yedg_orig[*,0] + skyoff
          wfslit.yedg_sky[*,1] = wfslit.yedg_orig[*,1] - skyoff
      endif
  endif

;  Check the Objects
  if not keyword_set( NOCHK ) then x_setobjgui, flux, wfobj, wfslit, wvimg=wave
  
; FITS
  ; Outfil
  objfil = 'Extract/Obj_'+strmid(wfccd[exp].img_root,3,3)+'.fits'
  wfccd[exp].obj_fil = objfil
  ; WRITE
  if not keyword_set(NOOBJ) then mwrfits, wfobj, objfil, /create
  if not keyword_set(NOSLIT) then mwrfits, wfslit, wfccd[exp].slit_fil, /create
   
  ; ALL DONE
  if not keyword_set(SILENT) then print, 'wfccd_mkobjstr: All Done!'
  return
end
