;+ 
; NAME:
; esi_lwdpltobj
;    Version 1.0
;
; PURPOSE:
;   Calls x_pltobj after gathering the relevant data
;
; CALLING SEQUENCE:
;   
;   esi_lwdpltobj, esi, objid, expsr, XSIZE=, YSIZE=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 1000)
;   YSIZE      - Size of gui in screen y-pixels (default = 600)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdpltobj, esi, objid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_lwdpltobj, esi, obj_id, obj_nm, expsr, XSIZE=xsize, $
                  YSIZE=ysize, FLUX=flux, FSPEC=fspec, XMAX=xmax, NOIMG=noimg

;
  if  N_params() LT 1 AND not keyword_set( FSPEC )  then begin 
    print,'Syntax - ' + $
      'esi_lwdpltobj, esi, obj_id, obj_nm, expsr, XSIZE=, YSIZE=, '
    print, '        /FLUX, /FSPEC, /NOIMG  [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then xsize = 1200
  if not keyword_set( XMAX ) and $
    (keyword_set(FSPEC) OR keyword_set(FLUX)) then xmax = 7.e-17

;;;;;;;
;  FSPEC

  if keyword_set( FSPEC ) then begin

      if N_params() EQ 0 then begin
          fils = findfile('Extract/*lwd.fits', count=nfil)
          if nfil EQ 0 then return 
          esi = x_guilist(fils)
      endif
      
      ; Read
      x_wrlwdfspec, lwdfspec, esi, /read

      ; Grab right object
      if keyword_set( obj_id) then begin
          indx = esi_lwdgetobjnm(lwdfspec, obj_id)
          objnm = obj_id
      endif else indx = x_getobjnm(lwdfspec, objnm)
      if indx EQ -1 then begin
          print, 'esi_lwdpltobj: Obj ', objnm, ' not found!'
          return
      endif

      ; Spectra
      spec_wv = lwdfspec[indx].wave[0:lwdfspec[indx].npix-1]
      spec_fx = lwdfspec[indx].fx[0:lwdfspec[indx].npix-1]
      spec_sig = fltarr(n_elements(spec_wv))
      a = where(lwdfspec[indx].var[0:lwdfspec[indx].npix-1] GT 0.)
      spec_sig[a] = float(sqrt(lwdfspec[indx].var[a]))


      ; Open Obj struct

      ; Open image file
      if not keyword_set( NOIMG ) then begin
          if not keyword_set( EXPSR ) then expsr = 0L
          esiobj = xmrdfits(lwdfspec[indx].obj_fil[expsr], 1, $
                          STRUCTYP='specobjstrct', /silent)
          obj = (where(esiobj.obj_id EQ lwdfspec[indx].obj_id, nobj))[0]
          if nobj NE 1 then begin
              print, 'wfccd_pltobj: Bad obj_id', nobj
              return
          endif
          wave = xmrdfits(esiobj[obj].slit_fil, /silent)
          fx = xmrdfits(esiobj[obj].spec2d_fil, 2, /silent)
      endif

      ; zabs
      if lwdfspec[indx].zans.z_err NE 0. then zin = lwdfspec[indx].zans.z
  endif else begin
; OBJ Alone ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Set Exposure
      if N_params() LT 2 then begin
          print,'Syntax - ' + $
            'esi_lwdpltobj, esi, obj_id, obj_nm, expsr, XSIZE=, '
          print, '        /FLUX, /FSPEC  [v1.0]'
          return
      endif
      allexp = where(esi.type EQ 'OBJ' AND esi.flg_anly NE 0 AND $
                     esi.obj_id EQ obj_id AND esi.mode EQ 1, nexp)
      if not keyword_set( EXPSR ) then expsr = 0L
      exp = allexp[expsr]

      ; Obj Structure
      esiobj = xmrdfits(esi[exp].obj_fil, 1, STRUCTYP='specobjstrct', /silent)


      ; Check for objnm
      if not keyword_set(obj_nm) then begin
          obj = x_getobjnm(esiobj, obj_nm)
      endif else begin
          obj = x_getobjnm(esiobj, obj_nm)
          objnm = strtrim(obj_nm,2)
      endelse

      ; Read fx, wave, var
      if not keyword_set( NOIMG ) then begin
          wave = xmrdfits(esi[exp].arc_fil, /silent)
;          var = mrdfits(esi[exp].img_final, 1, /silent)
          fx = xmrdfits(esi[exp].img_final, 2, /silent)
      endif

      ; Spectra
      npix = esiobj[obj].npix
      spec_wv = esiobj[obj].wave[0:npix-1]
      if keyword_set( FLUX ) then begin
          spec_fx = esiobj[obj].flux[0:npix-1] 
          spec_sig = esiobj[obj].sig[0:npix-1]
      endif else begin
          spec_fx = esiobj[obj].fx[0:npix-1]
          spec_sig = esiobj[obj].var[0:npix-1]
          a = where(esiobj[obj].var[0:npix-1] GT 0.)
          spec_sig[a] = sqrt(spec_sig[a])
      endelse

  endelse

  ;; Size
  if not keyword_set( NOIMG ) then begin
      sz = size(fx, /dimensions)

      ;; IMGWV
      imgwv = wave[lindgen(sz[0]), round(esiobj[obj].trace[lindgen(sz[0])])]
      mnwv = min(spec_wv, max=mxwv)
      mn = min(abs(mnwv-imgwv),imn)
      mn = min(abs(mxwv-imgwv),imx)
      imgwv = imgwv[imn:imx]

      ;; MASK ;;;;;;;;;;;;
      msk = bytarr(sz[0],sz[1])
      ymx = 0.
      ymn = 10000.
      for q=imn,imx do begin
          lmin = (esiobj[obj].trace[q] - 25) > 0L
          lmax = (esiobj[obj].trace[q] + 25) < (sz[1]-1)
          if lmax GT lmin then msk[q,lmin:lmax] = 1
          ymx = ymx > lmax
          ymn = ymn < lmin
      endfor
      badpix = where(msk EQ 0, nbad)
      if nbad NE 0 then begin
          fx[badpix] = 0.
          wave[badpix] = 0.
      endif
 
      ;; Get the sub images
      subfx = fx[imn:imx,ymn:ymx]
      subwv = wave[imn:imx,ymn:ymx]
      
      ;; x_pltobj
      x_pltobj, spec_wv, spec_fx, spec_sig, subfx, subwv, imgwv, $
        XSIZE=xsize, YSIZE=ysize, OBJNM=objnm, ZIN=zin, XMAX=xmax
  endif else begin
      x_specplot, spec_fx, spec_sig, WAVE=spec_wv, INFLG=4
  endelse

; GAME ON

  return
end
