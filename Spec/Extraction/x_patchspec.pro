;+ 
; NAME:
; x_combspec
;    Version 1.0
;
; PURPOSE:
;   Combines multiple exposures of the same slit
;
; CALLING SEQUENCE:
;   
;   x_combspec, wfccd, mask_id, exp_id
;
; INPUTS:
;   wfstrct     - WFCCD structure
;
; RETURNS:
;
; OUTPUTS:
;   wfspec      -  WFCCD fspec structure (fits file)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_combspec, wfccd, mask_id, exp_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_combspec_sngl, cnt, wffspec, wfobj, PATCH=patch

  wffspec[cnt].nexp = 1L
  wffspec[cnt].texp[0] = wfobj.exp

  ; Copy over
  tmpstr = wffspec[cnt]
  copy_struct, wfobj, tmpstr, EXCEPT=["wave", "fx", "var"]
  wffspec[cnt] = temporary(tmpstr)

  ; Add in obj_fil
  ipos = strpos(wfobj.spec2d_fil, 'ccd')
  obj_fil = 'Extract/Obj_'+strmid(wfobj.spec2d_fil, ipos+3)
  wffspec[cnt].obj_fil[0] = obj_fil

  ; Copy over flux, wave, var
  wffspec[cnt].wave = wfobj.wave[0:wfobj.npix-1]
  wffspec[cnt].fx = wfobj.flux[0:wfobj.npix-1]
  a = where(wfobj.sig GT 0., COMPLEMENT=badpix)
  wffspec[cnt].var[a] = double(wfobj.sig[a])^2
  
  ; Set wvmnx
  wffspec[cnt].wvmnx[0,0] = min(wfobj.wave[a], max=mx)
  wffspec[cnt].wvmnx[0,1] = mx

  ; Patch?
  if badpix[0] NE -1 then begin
      if keyword_set( PATCH ) then begin
          nbad = n_elements(badpix)
          for q=0L, nbad-1 do begin
              if badpix[q] GT wffspec[cnt].npix-1 then continue
              mn = 0 > (badpix[q]-5)
              mx = (wfobj.npix-1) < (badpix[q]+5)
              mdval = median(wfobj.flux[mn:mx])
              wffspec[cnt].fx[badpix[q]] = mdval
          endfor
      endif else wffspec[cnt].fx[badpix] = 0.
      wffspec[cnt].var[badpix] = -1.d
  endif

  return
end

;------------------------------------------------------------------------------
pro x_combspec_two, cnt, wffspec, wfobj, WVMNX=wvmnx, PATCH=patch

  if not keyword_set( WVMNX ) then wvmnx = [4500., 7000.]

  wffspec[cnt].nexp = 2L

  ; Check flux
  if wfobj[0].flg_flux NE wfobj[1].flg_flux then begin
      print, 'x_combspec_two: Wrong fluxes!'
      return
  endif

  ; Copy over
  tmpstr = wffspec[cnt]
  copy_struct, wfobj[0], tmpstr, EXCEPT=["wave", "fx", "var"]
  wffspec[cnt] = temporary(tmpstr)

  ; Exposure times
  wffspec[cnt].texp[0] = wfobj[0].exp
  wffspec[cnt].texp[1] = wfobj[1].exp

  ; wvmnx
  a = where(wfobj[0].sig GT 0.)
  wffspec[cnt].wvmnx[0,0] = min(wfobj[0].wave[a], max=mx)
  wffspec[cnt].wvmnx[0,1] = mx
  b = where(wfobj[1].sig GT 0.)
  wffspec[cnt].wvmnx[1,0] = min(wfobj[1].wave[b], max=mx)
  wffspec[cnt].wvmnx[1,1] = mx

  ; Add in obj_fil
  for q=0,1 do begin
      ipos = strpos(wfobj[q].spec2d_fil, 'ccd')
      obj_fil = 'Extract/Obj_'+strmid(wfobj[q].spec2d_fil, ipos+3)
      wffspec[cnt].obj_fil[q] = obj_fil
  endfor

  ; Set wavelength array
  wffspec[cnt].wave[0:wfobj[0].npix-1] = wfobj[0].wave[0:wfobj[0].npix-1]

  ; Normalize the flux
  gdwv = where(wfobj[0].wave GT wvmnx[0] AND $
               wfobj[0].wave LT wvmnx[1] AND $
               wfobj[0].sig GT 0. AND $
               wfobj[1].sig GT 0., ngd)
  if ngd EQ 0 then begin
      print, 'x_combspec_two: No way to normalize!'
      stop
      return
  endif
  rtio = wfobj[0].flux[gdwv] / wfobj[1].flux[gdwv]
  medflux = median(rtio)

  ; Reset the flux of second exposure
  
  fx2 = wfobj[1].flux * medflux
  sig2 = double(wfobj[1].sig * medflux)

;;;;;;;
; Add it up
               
  tmpfx = fltarr(wfobj[0].npix)
  tmpvar = replicate(-1.d, wfobj[0].npix) 

  ; Good pixels in both
  gdpix = where(wfobj[0].sig GT 0. AND wfobj[1].sig GT 0., ngd)
  tmpfx[gdpix] = (wfobj[0].flux[gdpix] + fx2[gdpix])/2.
  tmpvar[gdpix] = ( double(wfobj[0].sig[gdpix])^2 + sig2[gdpix]^2)/4.

  ; Good pixel in first
  onepix = where(wfobj[0].sig GT 0. AND wfobj[1].sig LE 0., ngd)
  if ngd GT 0 then begin
      tmpfx[onepix] = wfobj[0].flux[onepix] 
      tmpvar[onepix] = double(wfobj[0].sig[onepix])^2 
  endif

  ; Good pixel in the 2nd
  twopix = where(wfobj[1].sig GT 0. AND wfobj[0].sig LE 0., ngd)
  if ngd GT 0 then begin
      tmpfx[twopix] = wfobj[1].flux[twopix] 
      tmpvar[twopix] = sig2[twopix]^2 
  endif

  ; Good pixel in neither
  badpix = where(wfobj[0].sig LE 0. AND wfobj[1].sig LE 0., nbad)
  if nbad GT 0 then begin
      if keyword_set( PATCH ) then begin
          for q=0L, nbad-1 do begin
              if badpix[q] GT wffspec[cnt].npix-1 then continue
              mn = 0 > (badpix[q]-5)
              mx = (wffspec[cnt].npix-1) < (badpix[q]+5)
              mdval = median(wfobj[0].flux[mn:mx])
              tmpfx[badpix[q]] = mdval
          endfor
      endif else tmpfx[badpix] = 0.
  endif
      
;; FINISHING
  ; Set the flux
  wffspec[cnt].fx[0:wffspec[cnt].npix-1] = tmpfx
  wffspec[cnt].var[0:wffspec[cnt].npix-1] = tmpvar

  ; Set var values to -1
  a = where(wffspec[cnt].var[0:wffspec[cnt].npix-1] LE 0.)
  if a[0] NE -1 then wffspec[cnt].var[a] = -1.d

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

pro x_combspec, spec, var, PATCH=patch, OBJ_NM=OBJ_NM

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_combspec, wfccd, mask_id, [exp_id], /SILENT, OBJ_NM=, /PATCH [v1.0]'
    return
  endif 

;  Optional Keywords

; Set exp
  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp
  nexp = n_elements(exp)


;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'wfccd_combspec: Loading up the files...'

; Open Obj files
  for q=0L,nexp-1 do begin
      ; Read
      tmp = mrdfits(wfccd[exp[q]].obj_fil, 1, STRUCTYP='specobjstrct', /silent) 

      ; OBJ_NM
      if keyword_set(OBJ_NM) then begin
          gdobj = wfccd_getobjnm(tmp, obj_nm)
          if tmp[gdobj].flg_flux EQ 0 then begin
              print, 'wfccd_combspec: Data not good or not fluxed!'
              return
          endif 
      endif else begin
          ; Keep only good obj (require fluxed)
          gdobj = where(tmp.flg_flux NE 0 AND tmp.flg_anly NE 0, ngd)
          if ngd EQ 0 then begin
              print, 'wfccd_combspec: Data not good or not fluxed!'
              return
          endif
      endelse
      ; Add to total structure
      if q EQ 0 then wfobj = tmp[gdobj] else wfobj = [wfobj, tmp[gdobj]]
  endfor

  nspec = n_elements(wfobj)


;; ASCII OUTPUT ;;;

  if mask_id LT 10L then txtfil = 'Extract/Fspec_0'+$
    string(mask_id,format='(i1)')+'.txt' $
  else txtfil = 'Extract/Fspec_'+string(mask_id,format='(i2)')+'_comb.txt'
  close, /all
  if not keyword_set( OBJ_NM ) then openw, 56, txtfil

;;; CREATE FINAL ;;;

  if mask_id LT 10L then outfil = 'Extract/Fspec_0'+$
    string(mask_id,format='(i1)')+'.fits' $
  else outfil = 'Extract/Fspec_'+string(mask_id,format='(i2)')+'.fits'

  if not keyword_set( OBJ_NM ) then begin
      tmp = { wfccdfspecstrct }
      wffspec = replicate(tmp, nspec)
  endif else begin
      wfccd_wrfspec, wffspec, outfil, /read
  endelse

;;; FIND UNIQUE SLITS ;;;;

  if not keyword_set( OBJ_NM ) then begin
      uni_slit = uniq(wfobj.slit_id, sort(wfobj.slit_id))
      allslit = wfobj[uni_slit].slit_id
      nslit = n_elements(allslit)

  ; Loop on Slits

      cnt = 0L
      for q=0L, nslit-1 do begin
          slit = allslit[q]
          
      ; Find all unique obj
          allobj = where(wfobj.slit_id EQ slit)
          uni_obj = uniq(wfobj[allobj].obj_id, sort(wfobj[allobj].obj_id))
          objs = wfobj[allobj[uni_obj]].obj_id
          nobj = n_elements(objs)

      ; Loop on Obj
          for jj=0L, nobj-1 do begin
              obj = objs[jj]
              gdobj = where(wfobj.slit_id EQ slit AND wfobj.obj_id EQ obj $
                            AND wfobj.flg_anly NE 0, ngd)
              case ngd of 
                  1: wfccd_combspec_sngl, cnt, wffspec, wfobj[gdobj], PATCH=patch
                  2: wfccd_combspec_two, cnt, wffspec, wfobj[gdobj], PATCH=patch
                  else: stop
              endcase
          ; Output to ASCII
              printf, 56, wffspec[cnt].field+': '+$
                strtrim(wffspec[cnt].slit_id,2)+strtrim(wffspec[cnt].obj_id,2)
              for j=0L,wffspec[cnt].nexp-1 do begin
                  printf, 56, FORMAT='(10x,f7.1,1x,2f10.3,1x,a25)',$
                    wffspec[cnt].texp[j], $
                    wffspec[cnt].wvmnx[j,0], $
                    wffspec[cnt].wvmnx[j,1], $
                    wffspec[cnt].obj_fil[j]
              endfor
              cnt = cnt + 1
          endfor
      endfor
  endif else begin
      cnt = wfccd_getobjnm(wffspec, obj_nm)
      gdobj = lindgen(n_elements(wfobj))
      ngd = n_elements(wfobj)
      case ngd of 
          1: wfccd_combspec_sngl, cnt, wffspec, wfobj[gdobj], PATCH=patch
          2: wfccd_combspec_two, cnt, wffspec, wfobj[gdobj], PATCH=patch
          else: stop
      endcase
  endelse

;;;; OUTPUT  ;;;;

  if not keyword_set( OBJ_NM ) then wfccd_wrfspec, wffspec[0:cnt-1], outfil $
  else wfccd_wrfspec, wffspec, outfil

  close, /all

  print, 'wfccd_combspec:  All done!'


  return
end
  

