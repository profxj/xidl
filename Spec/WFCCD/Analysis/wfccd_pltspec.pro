;+ 
; NAME:
; wfccd_pltspec
;    Version 1.0
;
; PURPOSE:
;   Calls x_pltspec after gathering the relevant data
;
; CALLING SEQUENCE:
;   
;   wfccd_pltspec, wfccd, maskid, expsr, XSIZE=, YSIZE=
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
;   wfccd_pltspec, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_pltspec, wfccd, mask_id, expsr, obj_nm, XSIZE=xsize, $
                  YSIZE=ysize, FLUX=flux, FSPEC=fspec, XMAX=xmax

;
  if  N_params() LT 1 AND not keyword_set( FSPEC )  then begin 
    print,'Syntax - ' + $
      'wfccd_pltspec, wfccd, mask_id, expsr, obj_nm, XSIZE=, '
    print, '        /FLUX, /FSPEC  [v1.0]'
    return
  endif 

;  Optional Keywords

  if not keyword_set( XSIZE ) then xsize = 1200
  if not keyword_set( YSIZE ) then ysize = 800
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
      
      ; Read
      wfccd_wrfspec, wffspec, wfccd, /read

      ; Grab right object
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

      ; zabs
      if wffspec[indx].zans.z_err NE 0. then zin = wffspec[indx].zans.z
  endif else begin  ;; STANDARD
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

  if keyword_set( zin ) then begin
      synth = synthspec(wffspec[indx].zans, $
                        loglam=alog10(spec_wv))
  endif

  ;; Plot
  x_specplot, spec_fx, spec_sig, WAVE=spec_wv, /GAL, INFLG=4, $
    XSIZE=xsize, YSIZE=ysize, TITLE=objnm, ZIN=zin, $
    YTWO=synth


  return
end
