;+ 
; NAME:
; kast_pltobj
;    Version 1.1
;
; PURPOSE:
;   Calls x_pltobj after gathering the relevant data
;
; CALLING SEQUENCE:
;  kast_pltobj, kast, [setup, side, obj_id, obj_nm], EXPSR=, XSIZE=, $
;                 YSIZE=, /FLUX, /FSPEC, XMAX=, /NOIMG, /NOWV, /COMB
;
; INPUTS:
;   kast  --  Kast IDL structure
;  [setup]  --  Setup value
;   [side]  --  Specific camera [blue (1) vs. red (2)]
; [obj_id]  --  Object value
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOIMG  -- Do not display 2D image of spectrum
;  /COMB   -- Show the combined 1D spectra
;  EXPSR   -- Index of spectrum
;  /FLUX   -- Show the fluxed spectra
;  XMAX    -- Max height of spectrum [default: 7e-17]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   kast_pltobj, kast, /noimg
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Jul-2002 Written by JXP
;   ??  Major revisions by GEP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kast_pltobj, kast, setup, side, obj_id, obj_nm, EXPSR=expsr, XSIZE=xsize, $
                  YSIZE=ysize, FLUX=flux, FSPEC=fspec, XMAX=xmax, NOIMG=noimg, $
                 NOWV=nowv, COMB=comb, STD=std

;
  if  N_params() LT 1 AND not keyword_set( FSPEC )  then begin 
    print,'Syntax - ' + $
      'kast_pltobj, kast, setup, side, obj_id, [obj_nm], EXPSR=expsr, ' + $
      'XSIZE=, YSIZE=, '
    print, '        /FLUX, /FSPEC, /NOIMG  [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XMAX ) and $
    (keyword_set(FSPEC) OR keyword_set(FLUX)) then xmax = 7.e-17

;;;;;;;
;  FSPEC

  if keyword_set( FSPEC ) then begin

      if N_params() EQ 0 then begin
          fils = findfile('Extract/*.fits', count=nfil)
          if nfil EQ 0 then return 
          kast = x_guilist(fils)
      endif
      
      ; Read
      x_wrlwdfspec, lwdfspec, kast, /read

      ; Grab right object
      if keyword_set( obj_id) then begin
          indx = kast_lwdgetobjnm(lwdfspec, obj_id)
          objnm = obj_id
      endif else indx = x_getobjnm(lwdfspec, objnm)
      if indx EQ -1 then begin
          print, 'kast_lwdpltobj: Obj ', objnm, ' not found!'
          return
      endif

      ; Spectra
      if not keyword_set( NOWV ) then $
        spec_wv = lwdfspec[indx].wave[0:lwdfspec[indx].npix-1] $
      else spec_wv = dindgen(lwdfspec[indx].npix)
      spec_fx = lwdfspec[indx].fx[0:lwdfspec[indx].npix-1]
      spec_sig = fltarr(n_elements(spec_wv))
      a = where(lwdfspec[indx].var[0:lwdfspec[indx].npix-1] GT 0.)
      spec_sig[a] = float(sqrt(lwdfspec[indx].var[a]))


      ; Open Obj struct

      ; Open image file
      if not keyword_set( NOIMG ) then begin
          if not keyword_set( EXPSR ) then expsr = 0L
          kastobj = xmrdfits(lwdfspec[indx].obj_fil[expsr], 1, $
                          STRUCTYP='specobjstrct', /silent)
          obj = (where(kastobj.obj_id EQ lwdfspec[indx].obj_id, nobj))[0]
          if nobj NE 1 then begin
              print, 'wfccd_pltobj: Bad obj_id', nobj
              return
          endif
          wave = xmrdfits(kastobj[obj].slit_fil, /silent)
          fx = xmrdfits(kastobj[obj].spec2d_fil, 2, /silent)
      endif

      ; zabs
      if lwdfspec[indx].zans.z_err NE 0. then zin = lwdfspec[indx].zans.z
  endif else begin
; OBJ Alone ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ; Set Exposure
      if N_params() LT 4 then begin
          print,'Syntax - ' + $
            'kast_pltobj, kast, setup, side, obj_id, obj_nm, expsr, XSIZE=, '
          print, '        /FLUX, /FSPEC, /COMB  [v1.0]'
          return
      endif
      
      if not keyword_set(std) then begin
          allexp = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                         kast.side EQ side AND kast.setup EQ setup AND $
                         kast.obj_id EQ obj_id AND kast.type EQ 'OBJ', nexp)
          if not keyword_set( EXPSR ) then expsr = 0L
          exp = allexp[expsr]
      endif else begin
          allexp = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                         kast.side EQ side AND kast.setup EQ setup AND $
                         kast.type EQ 'STD', nexp)
          exp = allexp[obj_id]
      endelse 

      ; Obj Structure
      if not keyword_set(COMB) then begin 
        kastobj = xmrdfits(strtrim(kast[exp].obj_fil,2), 1, STRUCTYP='specobjstrct', /silent)

        ; Check for objnm
        if not keyword_set(obj_nm) then begin
            obj = x_getobjnm(kastobj, obj_nm)
        endif else begin
            obj = x_getobjnm(kastobj, obj_nm)
            objnm = strtrim(obj_nm,2)
        endelse

        ; Read fx, wave, var
        if not keyword_set( NOIMG ) then begin
;            var = mrdfits(kast[exp].img_final, 1, /silent)
            fx = xmrdfits(kast[exp].img_final, 0, /silent)
            sz = size(fx, /dimensions)
            wave = kastobj[obj].wave[0:kastobj[obj].npix-1] # replicate(1.,sz[1])
        endif

        ; Spectra
        npix = kastobj[obj].npix
        if not keyword_set( NOWV ) then spec_wv = kastobj[obj].wave[0:npix-1] $
        else spec_wv = dindgen(npix)
        if keyword_set( FLUX ) then begin
            spec_fx = kastobj[obj].flux[0:npix-1] 
            spec_sig = kastobj[obj].sig[0:npix-1]
        endif else begin
            spec_fx = kastobj[obj].fx[0:npix-1]
            spec_sig = kastobj[obj].var[0:npix-1]
            a = where(kastobj[obj].var[0:npix-1] GT 0.)
            spec_sig[a] = sqrt(spec_sig[a])
        endelse
      endif else begin
        ; it's a combined spectrum!
        infil = 'Extract/Fspec_'+kast[exp].img_root
;        kastobj = {kastspecstrct}
        kast_wrspec, kastobj, infil, /READ
        npix = kastobj.npix
        if not keyword_set( NOWV ) then spec_wv = kastobj.wave[0:npix-1] $
        else spec_wv = dindgen(npix)
        if keyword_set( FLUX ) then begin
            spec_fx = kastobj.flux[0:npix-1] 
            spec_sig = sqrt(kastobj.sig[0:npix-1])
        endif else begin
            spec_fx = kastobj.fx[0:npix-1]
            spec_sig = kastobj.var[0:npix-1]
            a = where(kastobj.var[0:npix-1] GT 0.)
            spec_sig[a] = sqrt(spec_sig[a])
        endelse
      endelse

  endelse

  ;; Size
  if not keyword_set( NOIMG ) AND not keyword_set(COMB) then begin
      sz = size(fx, /dimensions)

      ;; IMGWV
      imgwv = wave[lindgen(sz[0]), round(kastobj[obj].trace[lindgen(sz[0])])]
      mnwv = min(spec_wv, max=mxwv)
      mn = min(abs(mnwv-imgwv),imn)
      mn = min(abs(mxwv-imgwv),imx)
      imgwv = imgwv[imn:imx]

      ;; MASK ;;;;;;;;;;;;
      msk = bytarr(sz[0],sz[1])
      ymx = 0.
      ymn = 10000.
      for q=imn,imx do begin
          lmin = (kastobj[obj].trace[q] - 25) > 0L
          lmax = (kastobj[obj].trace[q] + 25) < (sz[1]-1)
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
