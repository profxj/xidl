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
                 STD=std, CHK=chk, FILE=file

;
  if  N_params() LT 2  AND NOT keyword_set( FILE) then begin 
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
      nobj = 1
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'esi_echflux: Loading up the files...'


; Open Flux Fit
  if not keyword_set( FXFIT ) then begin
;      fxfit = getenv('ESI_CALIBS')+'G191B2B_ECH.idl'
      fxfit = getenv('ESI_CALIBS')+'/G191B2B_ECH.idl'
      old = 1
  endif
  restore, fxfit

;;;;;;;
; LOOP+FLUX
  for qq=0L,nindx-1 do begin
      idx = indx[exp[qq]]
      ;; Open Obj file
      esiobj = xmrdfits(esi[indx[exp[qq]]].obj_fil, 1, /silent)

      if keyword_set( FORCE ) then begin
          esiobj.flg_flux = 0
          a = where(esiobj.flg_anly NE 0)
          esiobj[a].flg_anly = 1
      endif

      if keyword_set(NOBJ) then begin
          gdobj = lindgen(10*nobj)
      endif else begin
          gdobj = where(esiobj.order EQ 4 AND esiobj.flg_flux EQ 0, nobj)
      endelse
      
      for i=0L,nobj-1 do begin
          ;; Fluxing
          print, 'esi_echflux: Fluxing..'
          ;; Loop on Order
          obj_id = esiobj[gdobj[i]].obj_id
          for q=0L,9 do begin
              ;; Grab indx
              jj = where(esiobj.obj_id EQ obj_id AND esiobj.order EQ q)
              npix = esiobj[jj].npix
              if keyword_set( OLD ) then begin
                  newfx = x_fluxcalib(esiobj[jj].wave[0:npix-1], $
                                      esiobj[jj].fx[0:npix-1], flux_fit[q], $
                                      esiobj[jj].var[0:npix-1], newsig, $
                                      FLAMBDA=flambda)
                  ;; Take out exposure time
                  newsig = newsig / float(esi[indx[exp[qq]]].exp)
                  newfx = newfx / float(esi[indx[exp[qq]]].exp)
                  if keyword_set( CHK ) then $
                    x_splot, esiobj[jj].wave[0:npix-1], newfx, /block
                  esiobj[jj].flux[0:npix-1] = newfx
                  esiobj[jj].sig[0:npix-1] = newsig
              endif else begin
                  ;; New fluxing
                  kk = where(ordr_fit EQ (15-q))
                  gpx = lindgen(npix)
                  ;; Apply
                  full_rtio = x_calcfit(esiobj[jj].wave[gpx], $
                                        fitstr=tot_fit[kk])
                  newfx = esiobj[jj].fx[gpx]/full_rtio/ $
                    esi[idx].exp
;                  b = where(esiobj[jj].var[gpx] GT 0.)
                  newsig = sqrt(esiobj[jj].var[gpx]) $
                    / full_rtio / esi[idx].exp
                  esiobj[jj].flux[gpx] = newfx
                  esiobj[jj].sig[gpx] = newsig
              endelse
              ;; Save
              esiobj[jj].flg_anly = 2
              if keyword_set(FLAMBDA) then esiobj[jj].flg_flux = 2 $ ; flambda
              else esiobj[jj].flg_flux = 1 ; flambda
          endfor
      endfor
      ;; ;;
      ;; WRITE
      mwrfits, esiobj, esi[indx[exp[qq]]].obj_fil, /create
      spawn, 'gzip -f '+ esi[indx[exp[qq]]].obj_fil
  endfor


  print, 'esi_echflux: All Done!'

  return
end

