;+ 
; NAME:
; wfccd_skysub
;    Version 1.0
;
; PURPOSE:
;    Perform sky subtraction on the entire frame slit by slit
;
; CALLING SEQUENCE:
;   
;   wfccd_skysub, wfccd
;
; INPUTS:
;   wfccd     - WFCCD structure
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
;   wfccd_skysub, wfstrct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------
pro wfccd_skysub, wfccd, mask_id, exp_id, SLIT_ID=slit_id, NOFITS=nofits, $
                  DEBUG=debug, CHK=chk

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_skysub, wfccd, mask_id, [exp_id], SLIT_ID=, /NOFITS, /DEBUG' 
    print, '      /CHK  [v1.0]'
    return
  endif 

;  Optional Keywords

; Set exp
  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nexp)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp[0]

  if not keyword_set(SILENT) then $
    print, 'wfccd_skysub: Loading files...'

; Open Slit file
  slitstr = xmrdfits(wfccd[exp].slit_fil, 1, STRUCTYP='mslitstrct', /silent) 

; Open Obj file
  objstr = xmrdfits(wfccd[exp].obj_fil, 1, STRUCTYP='specobjstrct', /silent)

; Open Flux, Wave
  flux = xmrdfits(wfccd[exp].img_final, /silent)
  var  = xmrdfits(wfccd[exp].img_final, 1,/silent)
  wave = xmrdfits(wfccd[exp].img_final, 2, /silent)

;  Find good slits
  gdslit = where(slitstr.flg_anly NE 0, nslit)

; Set final images
  sz = size(flux, /dimensions)
  fin_img = fltarr(sz[0],sz[1])
;  sky_img = fltarr(sz[0],sz[1])
  all_rms = fltarr(sz[0], n_elements(slitstr))

; Open Prev sky as necessary
  if keyword_set(SLIT_ID) then begin
      fin_img = xmrdfits(wfccd[exp].img_final,3, /silent)
      all_rms = xmrdfits(wfccd[exp].img_final,4,/silent)
  endif
      
  if not keyword_set(SILENT) then $
    print, 'wfccd_skysub: Looping...'

;  LOOP 
  for q=0L,nslit-1 do begin
      ; Check for SLIT_ID
      if keyword_set(SLIT_ID) then begin
          if slitstr[gdslit[q]].id NE slit_id then continue
      endif
      if not keyword_set(SILENT) then print, 'Slit ', q, ' of ', nslit-1
      ; SUBTRACT
;      x_subskyslit, slitstr, gdslit[q], objstr, flux, wave, VAR=var, $
;        subimg=subimg, pix=pix, all_rms=rms, SKYIMG=subsky, DEBUG=debug
      wfccd_subskyslit, slitstr, gdslit[q], objstr, flux, wave, VAR=var, $
        subimg=subimg, pix=pix, all_rms=rms, SKYIMG=subsky, CHK=chk, $
        WVMNX=[3200., 11000.], DEBUG=debug
      fin_img[pix] = subimg[pix]
      all_rms[*,gdslit[q]] = rms
  endfor

; Output
  print, 'wfccd_skysub: Writing sky and RMS to '+wfccd[exp].img_final

  mwrfits, flux, wfccd[exp].img_final, /create
  mwrfits, var, wfccd[exp].img_final
  mwrfits, wave, wfccd[exp].img_final
  mwrfits, fin_img, wfccd[exp].img_final
  mwrfits, all_rms, wfccd[exp].img_final
  ;; COMPRESS
  print, 'wfccd_skysub: Compressing ', wfccd[exp].img_final
  spawn, 'gzip -f '+wfccd[exp].img_final

; All done
  print, 'wfccd_skysub: All done!'

  return
end
