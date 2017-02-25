;+ 
; NAME:
; wfccd_pltobj
;    Version 1.0
;
; PURPOSE:
;   Calls x_pltobj after gathering the relevant data
;
; CALLING SEQUENCE:
;   
;   wfccd_pltobj, wfccd, maskid, expsr, XSIZE=, YSIZE=
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
;   wfccd_pltobj, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_pltobj, wfccd, mask_id, expsr, obj_nm, XSIZE=xsize, $
                  YSIZE=ysize, FLUX=flux, FSPEC=fspec, XMAX=xmax

;
  if  N_params() LT 1 AND not keyword_set( FSPEC )  then begin 
    print,'Syntax - ' + $
      'wfccd_pltobj, wfccd, mask_id, expsr, obj_nm, XSIZE=, '
    print, '        /FLUX, /FSPEC  [v1.0]'
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
          fils = findfile('Extract/Fspec*fits', count=nfil)
          if nfil EQ 0 then return 
          wfccd = x_guilist(fils)
      endif
      
      ;; Read
      wfccd_wrfspec, wffspec, wfccd, /read

      ;; Grab right object
      if keyword_set( mask_id) then begin
          indx = x_getobjnm(wffspec, mask_id)
          objnm = mask_id
      endif else indx = x_getobjnm(wffspec, objnm)
      if indx EQ -1 then begin
          print, 'wfccd_pltobj: Obj ', objnm, ' not found!'
          return
      endif

      ; Spectra
      spec_wv = wffspec[indx].wave[0:wffspec[indx].npix-1]
      spec_fx = wffspec[indx].fx[0:wffspec[indx].npix-1]
      spec_sig = fltarr(n_elements(spec_wv))
      a = where(wffspec[indx].var[0:wffspec[indx].npix-1] GT 0.)
      spec_sig[a] = float(sqrt(wffspec[indx].var[a]))

      ; Open Obj struct
      if not keyword_set( EXPSR ) then expsr = 0L
      wfobj = xmrdfits(strtrim(wffspec[indx].obj_fil[expsr],2), 1, $
                      STRUCTYP='specobjstrct', /silent)
      obj = (where(wffspec[indx].slit_id EQ wfobj.slit_id AND $
                   wfobj.obj_id EQ wffspec[indx].obj_id, nobj))[0]
      if nobj NE 1 then begin
          print, 'wfccd_pltobj: Bad obj_id', nobj
          return
      endif

      ; Open image file
      wave = xmrdfits(wfobj[obj].spec2d_fil, 2, /silent)
      fx = xmrdfits(wfobj[obj].spec2d_fil, 3, /silent)

      ; zabs
      if wffspec[indx].zans.z_err NE 0. then zin = wffspec[indx].zans.z
  endif else begin
      ; Set Exposure
      if N_params() LE 2 then begin
          print,'Syntax - ' + $
            'wfccd_pltobj, wfccd, mask_id, expsr, obj_nm, XSIZE=, '
          print, '        /FLUX, /FSPEC  [v1.0]'
          return
      endif
      allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
                     wfccd.mask_id EQ mask_id, nexp)
      exp = allexp[expsr]

      ; Obj Structure
      wfobj = xmrdfits(wfccd[exp].obj_fil, 1, STRUCTYP='specobjstrct', /silent)

      ; Check for objnm
      if not keyword_set(obj_nm) then begin
          obj = x_getobjnm(wfobj, objnm)
      endif else begin
          obj = x_getobjnm(wfobj, obj_nm)
          objnm = strtrim(obj_nm,2)
      endelse

      ; Read fx, wave, var
      wave = xmrdfits(wfccd[exp].img_final, 2, /silent)
      var = xmrdfits(wfccd[exp].img_final, 1, /silent)
      fx = xmrdfits(wfccd[exp].img_final, 3, /silent)

      ; Spectra
      spec_wv = wfobj[obj].wave
      if keyword_set( FLUX ) then begin
          spec_fx = wfobj[obj].flux 
          spec_sig = wfobj[obj].sig
      endif else begin
          spec_fx = wfobj[obj].fx
          spec_sig = wfobj[obj].var
          a = where(wfobj[obj].var GT 0.)
          spec_sig[a] = sqrt(spec_sig[a])
      endelse

  endelse

  ;;;;;;;; SLIT ;;;;;;;;;;;;;
  slitstr = xmrdfits(wfobj[obj].slit_fil, 1, /silent)
;  slitstr = xmrdfits(wfobj[obj].slit_fil, 1, STRUCTYP='mslitstrct', /silent)
  slit = where(slitstr.id EQ wfobj[obj].slit_id) 

  sz = size(fx, /dimensions)

  ;; KLUDGE
  if max(slitstr[slit].yedg_sky[*,0]) LT 1. then begin
      slitstr[slit].yedg_sky[*,0] = slitstr[slit].yedg_orig[*,0]
      slitstr[slit].yedg_sky[*,1] = slitstr[slit].yedg_orig[*,1]
  endif

  ;;;;;;;; MASK ;;;;;;;;;;;;
  msk = bytarr(sz[0],sz[1])
  ymx = 0.
  ymn = 10000.
  for q=0L,sz[0]-1 do begin
      lmin = (wfobj[obj].trace[q] - 25) > round(slitstr[slit].yedg_sky[q,0])
      lmax = (wfobj[obj].trace[q] + 25) < round(slitstr[slit].yedg_sky[q,1])
      if lmax GT lmin then msk[q,lmin:lmax] = 1
      ymx = ymx > lmax
      ymn = ymn < lmin
  endfor
  badpix = where(msk EQ 0, nbad)
  if nbad NE 0 then begin
      fx[badpix] = 0.
      wave[badpix] = 0.
  endif
 
  ; Get the sub images
  subfx = fx[*,ymn:ymx]
  subwv = wave[*,ymn:ymx]
;  subvar = var[*,ymn:ymx]

;;;;;;;;;;;;;
; FLIP
  subfx = rotate(subfx, 5) ; Flip x
  subwv = rotate(subwv, 5) ; Flip x
;  subvar = rotate(subwv, 5) ; Flip x

;;;;;;;;;;;;;
; IMG WAVE
  dumwv = wave[lindgen(sz[0]), round(wfobj[obj].trace[lindgen(sz[0])])]
  imgwv = rotate(dumwv,2)
  delvarx, dumwv

; GAME ON

  x_pltobj, spec_wv, spec_fx, spec_sig, subfx, subwv, imgwv, $
    XSIZE=xsize, YSIZE=ysize, OBJNM=objnm, ZIN=zin, XMAX=xmax

  return
end
