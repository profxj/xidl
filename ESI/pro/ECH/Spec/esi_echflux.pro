;+ 
; NAME:
; esi_echflux
;    Version 1.0
;
; PURPOSE:
;   Fluxes esi data.  Puts in flambda units by default
;
; CALLING SEQUENCE:
;   
;   esi_echflux, esi, obj_id, [exp_id]
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
;   esi_echflux, esi, obj_id, exp_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro esi_echflux, esi, obj_id, exp_id, FXFIT=fxfit, FNU=fnu, FORCE=force, $
                 STD=std, CHK=chk

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'esi_echflux, esi, obj_id, [exp_id], /FNU, /FORCE, FXFIT=, /STD [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( FNU ) then flambda = 1

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                   esi.obj_id EQ obj_id AND $
                   strtrim(esi.type,2) EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'esi_echskysub: No images to sky subtract!', obj_id
          return
      endif
  endif else begin
      indx = obj_id[0]
      nindx = 1L
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'esi_echflux: Loading up the files...'


; Open Flux Fit
  if not keyword_set( FXFIT ) then begin
      fxfit = getenv('XIDL_DIR')+'/ESI/CALIBS/G191B2B_ECH.idl'
      restore, fxfit
  endif

;;;;;;;
; LOOP+FLUX
  for qq=0L,nindx-1 do begin
      ;; Open Obj file
      esiobj = mrdfits(esi[indx[exp[qq]]].obj_fil, 1, $
                       STRUCTYP='dblsobjstrct', /silent)

      if keyword_set( FORCE ) then begin
          esiobj.flg_flux = 0
          a = where(esiobj.flg_anly NE 0)
          esiobj[a].flg_anly = 1
      endif
      gdobj = where(esiobj.slit_id EQ 4 AND esiobj.flg_flux EQ 0, nobj)

      for i=0L,nobj-1 do begin
          ;; Loop on Order
          obj_id = esiobj[gdobj[i]].obj_id
          for q=0L,9 do begin
              ;; Grab indx
              jj = where(esiobj.obj_id EQ obj_id AND esiobj.slit_id EQ q)
              npix = esiobj[jj].npix
              newfx = x_fluxcalib(esiobj[jj].wave[0:npix-1], $
                                  esiobj[jj].fx[0:npix-1], flux_fit[q], $
                                  esiobj[jj].var[0:npix-1], newsig, $
                                  FLAMBDA=flambda)
              ;; Take out exposure time
              newsig = newsig / float(esi[indx[exp[qq]]].exp)
              newfx = newfx / float(esi[indx[exp[qq]]].exp)
              if keyword_set( CHK ) then $
                x_splot, esiobj[jj].wave[0:npix-1], newfx, /block
              ;; Save
              esiobj[jj].flg_anly = 2
              if keyword_set(FLAMBDA) then esiobj[jj].flg_flux = 2 $ ; flambda
              else esiobj[jj].flg_flux = 1 ; flambda
              esiobj[jj].flux[0:npix-1] = newfx
              esiobj[jj].sig[0:npix-1] = newsig
          endfor
      endfor
      ;; ;;
      ;; WRITE
      mwrfits, esiobj, esi[indx[exp[qq]]].obj_fil, /create
  endfor


  print, 'esi_echflux: All Done!'

  return
end

