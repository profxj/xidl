;+ 
; NAME:
; esi_lwdflux
;    Version 1.0
;
; PURPOSE:
;   Fluxes esi data.  Puts in flambda units by default
;
; CALLING SEQUENCE:
;   
;   esi_lwdflux, esi, obj_id, [exp_id]
;
; INPUTS:
;   esi     - WFCCD structure
;   obj_id     - WFCCD structure
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
;   esi_lwdflux, esi, obj_id, exp_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro esi_lwdflux, esi, obj_id, exp_id, FXFIT=fxfit, FNU=fnu, FORCE=force

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'esi_lwdflux, esi, obj_id, [exp_id], /FNU, /FORCE, FXFIT= [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( FNU ) then flambda = 1

; Set exp
  allexp = where(esi.type EQ 'OBJ' AND esi.flg_anly NE 0 AND $
              esi.obj_id EQ obj_id AND esi.mode EQ 1)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp
  nexp = n_elements(exp)

;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'esi_lwdflux: Loading up the files...'


; Open Flux Fit
  if not keyword_set( FXFIT ) then $
    fxfit = getenv('XIDL_DIR')+'/ESI/pro/LWD/Stds/bd284211_fit.fits'
  x_fitstrtofits, flux_fit, fxfit, /reverse

;;;;;;;
; LOOP+FLUX
  for qq=0L,nexp-1 do begin
      esiobj = xmrdfits(esi[exp[qq]].obj_fil, 1, STRUCTYP='specobjstrct', /silent)

      if keyword_set( FORCE ) then begin
          esiobj.flg_flux = 0
          a = where(esiobj.flg_anly NE 0)
          esiobj[a].flg_anly = 1
      endif
      gdobj = where(esiobj.flg_anly NE 0 AND esiobj.flg_flux EQ 0, nobj)

      for i=0L,nobj-1 do begin
          ;; Open Obj file
          npix = esiobj[gdobj[i]].npix
          newfx = x_fluxcalib(esiobj[gdobj[i]].wave[0:npix-1], $
                              esiobj[gdobj[i]].fx[0:npix-1], flux_fit, $
                              esiobj[gdobj[i]].var[0:npix-1], newsig, $
                              FLAMBDA=flambda)
          ;; Take out exposure time
          newsig = newsig / float(esi[exp[qq]].exp)
          newfx = newfx / float(esi[exp[qq]].exp)
          ;; Save
          esiobj[gdobj[i]].flg_anly = 2
          if keyword_set(FLAMBDA) then esiobj[gdobj[i]].flg_flux = 2 $; flambda
          else esiobj[gdobj[i]].flg_flux = 1 ; flambda
          esiobj[gdobj[i]].flux[0:npix-1] = newfx
          esiobj[gdobj[i]].sig[0:npix-1] = newsig
      endfor
      ;; ;;
      ;; WRITE
      mwrfits, esiobj, esi[exp[qq]].obj_fil, /create
  endfor


  print, 'esi_lwdflux: All Done!'

  return
end

