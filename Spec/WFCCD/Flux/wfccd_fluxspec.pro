;+ 
; NAME:
; wfccd_fluxspec
;    Version 1.0
;
; PURPOSE:
;   Fluxes wfccd data.  Puts in flambda units by default
;
; CALLING SEQUENCE:
;   
;   wfccd_fluxspec, wfccd, mask_id, [exp_id]
;
; INPUTS:
;   wfccd     - WFCCD structure
;   mask_id     - WFCCD structure
;   [exp_id]     - WFCCD structure
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
;   wfccd_fluxspec, wfccd, mask_id, exp_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro wfccd_fluxspec, wfccd, mask_id, exp_id, FXFIT=fxfit, FNU=fnu, FORCE=force

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_fluxspec, wfccd, mask_id, [exp_id], FNU=fnu [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( FNU ) then FLAMBDA=1

; Set exp
  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nexp)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp[0]

;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'wfccd_fluxspec: Loading up the files...'

; Open Slit file
  wfslit = xmrdfits(wfccd[exp].slit_fil, 1, STRUCTYP='mslitstrct', /silent) 

; Open Obj file
  wfobj = xmrdfits(wfccd[exp].obj_fil, 1, STRUCTYP='specobjstrct', /silent)

; Open Flux Fit
  if not keyword_set( FXFIT ) then $
    fxfit = getenv('XIDL_DIR')+'/Spec/WFCCD/Flux/Stds/LTT7379_fit.fits'
  x_fitstrtofits, flux_fit, fxfit, /reverse

;;;;;;;
; LOOP+FLUX

  if keyword_set( FORCE ) then begin
      wfobj.flg_flux = 0
      a = where(wfobj.flg_anly NE 0)
      wfobj[a].flg_anly = 1
  endif
  gdobj = where(wfobj.flg_anly EQ 1 AND wfobj.flg_flux EQ 0, nobj)

  for i=0L,nobj-1 do begin
      npix = wfobj[gdobj[i]].npix
      newfx = x_fluxcalib(wfobj[gdobj[i]].wave[0:npix-1], $
                          wfobj[gdobj[i]].fx[0:npix-1], flux_fit, $
                          wfobj[gdobj[i]].var[0:npix-1], newsig, $
                          FLAMBDA=flambda, $
                          TRUCONV=truconv)
      ; Take out exposure time
      newsig = newsig / float(wfobj[gdobj[i]].exp)
      newfx = newfx / float(wfobj[gdobj[i]].exp)

      ; Save
      wfobj[gdobj[i]].flg_anly = 2
      wfobj[gdobj[i]].flg_flux = 2  ; flambda
      wfobj[gdobj[i]].flux[0:npix-1] = newfx
      wfobj[gdobj[i]].sig[0:npix-1] = newsig
  endfor

;;;;
; WRITE
  mwrfits, wfobj, wfccd[exp].obj_fil, /create

  print, 'wfccd_fluxspec: All Done!'

  return
end

