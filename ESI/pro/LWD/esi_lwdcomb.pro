;+ 
; NAME:
; esi_lwdcomb
;    Version 1.0
;
; PURPOSE:
;   Combines multiple spectral exposures of the same obj
;
; CALLING SEQUENCE:
;   
;   esi_lwdcomb, esi, obj_id, exp_id
;
; INPUTS:
;   esi     - WFCCD structure
;
; RETURNS:
;
; OUTPUTS:
;   lwdfspec      -  WFCCD fspec structure (fits file)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdcomb, esi, obj_id, exp_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro esi_lwdcomb_out, lwdfspec

  ;; Output to ASCII
  printf, 56, lwdfspec.field
  for j=0L,lwdfspec.nexp-1 do begin
      printf, 56, FORMAT='(10x,f7.1,1x,2f10.3,1x,a25)',$
        lwdfspec.texp[j], $
        lwdfspec.wvmnx[j,0], $
        lwdfspec.wvmnx[j,1], $
        lwdfspec.obj_fil[j]
  endfor
  return
end

;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro esi_lwdcomb, esi, obj_id, exp_id, SILENT=silent, PATCH=patch, $
                    OBJ_NM=OBJ_NM

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'esi_lwdcomb, esi, obj_id, [exp_id], /SILENT, OBJ_NM=, /PATCH [v1.0]'
    return
  endif 

;  Optional Keywords

; Set exp
  allexp = where(esi.type EQ 'OBJ' AND esi.flg_anly NE 0 AND $
              esi.mode EQ 1 AND esi.obj_id EQ obj_id)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp
  nexp = n_elements(exp)


;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'esi_lwdcomb: Loading up the files...'

; Open Obj files
  for q=0L,nexp-1 do begin
      ; Read
      tmp = xmrdfits(esi[exp[q]].obj_fil, 1, STRUCTYP='specobjstrct', /silent) 

      ; OBJ_NM
      if keyword_set(OBJ_NM) then begin
          gdobj = esi_getobjnm(tmp, obj_nm)
          if tmp[gdobj].flg_flux EQ 0 then begin
              print, 'esi_lwdcomb: Data not good or not fluxed!'
              return
          endif 
      endif else begin
          ; Keep only good obj (require fluxed)
          gdobj = where(tmp.flg_flux NE 0 AND tmp.flg_anly NE 0, ngd)
          if ngd EQ 0 then begin
              print, 'esi_lwdcomb: Data not good or not fluxed!'
              return
          endif
      endelse
      ; Add to total structure
      if q EQ 0 then esiobj = tmp[gdobj] else esiobj = [esiobj, tmp[gdobj]]
  endfor

  nspec = n_elements(esiobj)


;; ASCII OUTPUT ;;;

  if obj_id LT 10L then txtfil = 'Extract/LWDF_0'+$
    string(obj_id,format='(i1)')+'.txt' $
  else txtfil = 'Extract/LWDF_'+string(obj_id,format='(i2)')+'.txt'
  close, /all
  if not keyword_set( OBJ_NM ) then openw, 56, txtfil

;;; CREATE FINAL ;;;

  outfil = 'Extract/'+strtrim(esi[exp[0]].Obj,2)+'_lwd.fits'

  if not keyword_set( OBJ_NM ) then begin
      tmp = { lwdfspecstrct }
      lwdfspec = replicate(tmp, nspec)
  endif else begin
      esi_wrfspec, lwdfspec, outfil, /read
  endelse

;;; LOOP ON OBJ;;;;

  if not keyword_set( OBJ_NM ) then begin

      cnt = 0L
      ;; Science first
      sciobj = where(esiobj.obj_id EQ 'a' AND esiobj.flg_anly NE 0, nsci)
      if nsci NE 0 then begin
          ;; Copy
          lwdfspec[cnt].nexp = nsci
          for i=0L,nsci-1 do lwdfspec[cnt].texp[i] = esiobj[sciobj[i]].exp
          tmpstr = lwdfspec[cnt]
          copy_struct, esiobj[sciobj[0]], tmpstr, EXCEPT=["wave","fx","var"]
          lwdfspec[cnt] = tmpstr
          for i=0L,nsci-1 do begin
              ipos = strpos(esiobj[sciobj[i]].spec2d_fil, 'esi')
              obj_fil = 'Extract/Obj_'+strmid(esiobj[sciobj[i]].spec2d_fil, $
                                              ipos)
              lwdfspec[cnt].obj_fil[i] = obj_fil
          endfor
          ;; NPIX
          npix = lwdfspec[cnt].npix
          ;; Coadd
          if nsci EQ 1 then begin
              lwdfspec[cnt].wave[0:npix-1] = esiobj[sciobj].wave[0:npix-1]
              lwdfspec[cnt].fx[0:npix-1] = esiobj[sciobj].flux[0:npix-1]
              lwdfspec[cnt].var[0:npix-1] = $
                double(esiobj[sciobj].sig[0:npix-1])^2
          endif else begin
              var = double(esiobj[sciobj].sig[0:npix-1])^2
              for kk=0L,nsci-1 do begin
                  a = where(esiobj[sciobj[kk]].sig[0:npix-1] LT 0., na)
                  if na GT 0 then var[a,kk] = -1.
              endfor
              x_combspec, esiobj[sciobj].flux[0:npix-1], var, $
                fflux, fvar, WAVE=esiobj[sciobj[0]].wave[0:npix-1], $
                NRMFLUX=[4500., 7000.]
              lwdfspec[cnt].wave[0:npix-1] = esiobj[sciobj[0]].wave[0:npix-1]
              lwdfspec[cnt].fx[0:npix-1] = temporary(fflux[0:npix-1])
              lwdfspec[cnt].var[0:npix-1] = temporary(fvar[0:npix-1])
          endelse
          ;; Output
          esi_lwdcomb_out, lwdfspec[cnt]
          ;;
          cnt = cnt+1
      endif
                  
      ;; SERENDIP
      sdpobj = where(esiobj.obj_id NE 'a' AND esiobj.flg_anly NE 0, nsdp)
      if nsdp NE 0 then begin
          srt = sort(esiobj[sdpobj].ycen)
          jj = 0L
          kk= 0L
          while(jj LT nsdp) do begin
              ;; Find all within 2pix
              kp = where( abs(esiobj[sdpobj].ycen-esiobj[sdpobj[srt[jj]]].ycen) $
                          LT 2., ngd)
              jj = jj+ngd
              kk = kk+1
              ;; Reset
              gd = sdpobj[kp]
              ;; Copy
              lwdfspec[cnt].nexp = ngd
              for i=0L,ngd-1 do lwdfspec[cnt].texp[i] = esiobj[gd[i]].exp
              tmpstr = lwdfspec[cnt]
              copy_struct, esiobj[gd[0]], tmpstr, EXCEPT=["wave","fx","var"]
              lwdfspec[cnt] = tmpstr
              for i=0L,ngd-1 do begin
                  ipos = strpos(esiobj[gd[i]].spec2d_fil, 'esi')
                  obj_fil = 'Extract/Obj_'+strmid(esiobj[gd[i]].spec2d_fil, ipos+3)
                  lwdfspec[cnt].obj_fil[i] = obj_fil
              endfor
              ;; OBJ name
              lwdfspec[cnt].obj_id = x_objnumid(kk)
              ;; NPIX
              npix = lwdfspec[cnt].npix
              ;; Coadd
              if ngd EQ 1 then begin
                  lwdfspec[cnt].wave[0:npix-1] = esiobj[gd].wave[0:npix-1]
                  lwdfspec[cnt].fx[0:npix-1] = esiobj[gd].flux[0:npix-1]
                  gdpix = where(esiobj[gd].sig[0:npix-1] GT 0., COMPLEMENT=badpix, $
                                NCOMPLEMENT=nbad)
                  lwdfspec[cnt].var[gdpix] = double(esiobj[gd].sig[gdpix])^2
                  if nbad NE 0 then $
                    lwdfspec[cnt].var[badpix] = double(esiobj[gd].sig[badpix])
              endif else begin
                  var = double(esiobj[gd].sig[0:npix-1])^2
                  for kk=0L,ngd-1 do begin
                      a = where(esiobj[gd[kk]].sig[0:npix-1] LT 0., na)
                      if na GT 0 then var[a,kk] = -1.
                  endfor
                  x_combspec, esiobj[gd].flux[0:npix-1], var, $
                    fflux, fvar, WAVE=esiobj[gd[0]].wave[0:npix-1], $
                    NRMFLUX=[5000., 7000.]
                  lwdfspec[cnt].wave[0:npix-1] = esiobj[gd[0]].wave[0:npix-1]
                  lwdfspec[cnt].fx[0:npix-1] = temporary(fflux[0:npix-1])
                  lwdfspec[cnt].var[0:npix-1] = temporary(fvar[0:npix-1])
              endelse
              cnt = cnt+1
              ;; Output
              esi_lwdcomb_out, lwdfspec[cnt]
          endwhile
      endif
  endif else begin         ;;;;;; OBJ_NM ;;;;;;
      slen = strlen(obj_nm)
      obj = strmid(obj_nm,slen-1)
      gd = where(esiobj.obj_id EQ obj AND esiobj.flg_anly NE 0, ngd)
      if ngd EQ 0 then begin
          print, 'esi_lwdcomb: No spectra found!'
          return
      endif else begin
          ;; Set cnt
          cnt = where(lwdfspec.obj_id EQ obj )
          ;; Copy
          lwdfspec[cnt].nexp = ngd
          for i=0L,ngd-1 do lwdfspec[cnt].texp[i] = esiobj[gd[i]].exp
          tmpstr = lwdfspec[cnt]
          copy_struct, esiobj[gd[0]], tmpstr, EXCEPT=["wave","fx","var"]
          lwdfspec[cnt] = tmpstr
          for i=0L,ngd-1 do begin
              ipos = strpos(esiobj[gd[i]].spec2d_fil, 'esi')
              obj_fil = 'Extract/Obj_'+strmid(esiobj[gd[i]].spec2d_fil, ipos+3)
              lwdfspec[cnt].obj_fil[i] = obj_fil
          endfor
          ;; NPIX
          npix = lwdfspec[cnt].npix
          ;; Coadd
          if ngd EQ 1 then begin
              lwdfspec[cnt].wave[0:npix-1] = esiobj[gd].wave[0:npix-1]
              lwdfspec[cnt].fx[0:npix-1] = esiobj[gd].flux[0:npix-1]
              gdpix = where(esiobj[gd].sig[0:npix-1] GT 0., COMPLEMENT=badpix, $
                            NCOMPLEMENT=nbad)
              lwdfspec[cnt].var[gdpix] = double(esiobj[gd].sig[gdpix])^2
              if nbad NE 0 then $
                lwdfspec[cnt].var[badpix] = double(esiobj[gd].sig[badpix])
          endif else begin
              var = double(esiobj[gd].sig[0:npix-1])^2
              for kk=0L,ngd-1 do begin
                  a = where(esiobj[gd[kk]].sig[0:npix-1] LT 0., na)
                  if na NE 0 then var[a,kk] = -1.
              endfor
              x_combspec, esiobj[gd].flux[0:npix-1], var, $
                fflux, fvar, WAVE=esiobj[gd[0]].wave[0:npix-1], $
                NRMFLUX=[4500., 7000.]
              lwdfspec[cnt].wave[0:npix-1] = esiobj[gd[0]].wave[0:npix-1]
              lwdfspec[cnt].fx[0:npix-1] = temporary(fflux)
              lwdfspec[cnt].var[0:npix-1] = temporary(fvar)
          endelse
          cnt = cnt+1
          ;; Output
;          esi_lwdcomb_out, lwdfspec[cnt]
      endelse
  endelse
              

;;;; OUTPUT  ;;;;

  if not keyword_set( OBJ_NM ) then x_wrlwdfspec, lwdfspec[0:cnt-1], outfil $
  else x_wrlwdfspec, lwdfspec, outfil

  close, /all

  print, 'esi_lwdcomb:  All done!'


  return
end
  

