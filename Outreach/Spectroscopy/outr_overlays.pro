;+ 
; NAME:
; outr_overlays
;   Version 1.1
;
; PURPOSE:
;    GUI for plotting a spectrum and doing quick analysis 
;
; CALLING SEQUENCE:
;   
;   outr_spectratool_ flux, [ysin], XSIZE=, YSIZE=, TITLE=, WAVE=, LLIST=,
;           /QAL, /GAL, ZIN=, /BLOCK, /NRM, /LLS
;
; INPUTS:
;   flux  - Flux array (or FITS file)
;   [ysin]  - Sigma array (or FITS file)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave=      - wavelength array
;   dunit=     - Extension in the FITS file
;   INFLG=     - Specifies the type of input files (or arrays)
;   /LLS       - Use Lyman limit line list
;   /QSO       - Use Quasar line list
;   /GAL       - Use galaxy line list
;   /QAL       - Use quasar absorption line list
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   outr_spectratool_ 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   01-Apr-2008 Written by JXP
;-
;------------------------------------------------------------------------------

function outr_overlays_mkwv, wvmnx, npix
  if wvmnx[1] GT wvmnx[0]*10 then begin ; Log
      wv = 10.^(alog10(wvmnx[0]) + findgen(npix)*$
                           (alog10(wvmnx[1]/wvmnx[0]))/float(npix-1))
  endif else begin
      wv = wvmnx[0] + findgen(npix)*(wvmnx[1]-wvmnx[0])/float(npix-1)
  endelse

return, wv
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro outr_overlays_lamp, strct, indx

  npix = n_elements(strct.wave)
  strct.npix = npix

  ;; Lamps
  case indx of
      0: begin  ; Na Lamp
          sub_wv = 5890. + findgen(npix-2)*10/float(npix-2)
          strct.wave = [1e-10, sub_wv, 1e10]
          sig = 5895.*5./3e5 ; Ang
          sub_fx = exp(-1*abs(sub_wv-5891.5833d)/(2*sig^2))*2 + $
                   exp(-1*abs(sub_wv-5897.5581d)/(2*sig^2))
          mx = max(sub_fx)
          strct.fx = [0., sub_fx/mx, 0.]
          strct.f_norm = 1.
      end
      else: stop
  endcase
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro outr_overlays_star, strct, INIT=init, WVMNX=iwvmnx, PARAM=param, $
  FIDDLE=fiddle

  npix = n_elements(strct.wave)
  strct.npix = npix

  if keyword_set(IWVMNX) then wvmnx = iwvmnx
  if keyword_set(PARAM) then strct.param[0:n_elements(param)-1] = param

  ;; Initialize
  if keyword_set(INIT) then begin
      if not keyword_set(PARAM) then strct.param[0] = 5500.  ; Temperature
      ;; Use Wien's Law for wavelength range (if not specified)
      if not keyword_set(WVMNX) then begin
          wien = 2.9e7 / strct.param[0]  ;; Ang
          wvmnx = [wien/10., wien*10]
      endif
      strct.wave = outr_overlays_mkwv(wvmnx, npix)
      dblf = blackbody(strct.wave, strct.param[0]) 
      strct.f_norm = max(dblf)
      strct.fx = float(dblf/strct.f_norm)
;      if strmatch(strct.name,'Human') then stop
  endif

  if keyword_set(FIDDLE) then begin
      if keyword_set(WVMNX) then $
        strct.wave = outr_overlays_mkwv(wvmnx, npix)
      dblf = blackbody(strct.wave, strct.param[0]) 
      strct.f_norm = max(dblf)
      strct.fx = float(dblf/strct.f_norm)
  endif

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro outr_overlays_qso, strct, INIT=init, PARAM=param, $
  FIDDLE=fiddle

  npix = n_elements(strct.wave)
  strct.npix = npix

  if keyword_set(IWVMNX) then wvmnx = iwvmnx
  if keyword_set(PARAM) then strct.param[0:n_elements(param)-1] = param

  ;; Initialize
  if keyword_set(INIT) then begin
      if not keyword_set(PARAM) then strct.param[0] = 5500.  ; Temperature
      ;; Read in spectrum
      fx = x_readspec(getenv('XIDL_DIR')+ $
                      '/Outreach/Spectroscopy/Spectra/vanden.fits', wav=wav, $
                      inflg=2, npix=n2)
      strct.wave[0:n2-1] = wav
      strct.f_norm = max(fx)
      strct.fx[0:n2-1] = float(fx/strct.f_norm)
      strct.wave[n2:*] = strct.wave[n2-1]
;      if strmatch(strct.name,'Human') then stop
  endif

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro outr_overlays, strct, INIT=init, FIDDLE=fiddle, PARAM=param, $
                   WVMNX=wvmnx

;
  if  N_params() LT 1 and not keyword_set(DEFAULT) then begin 
    print,'Syntax - ' + $
             'outr_overlay, strct, ) [v1.0]'
    return
  endif 

  gd = where(strlen(strct.name) GT 0, nstrct)
  if nstrct EQ 0 then return

  ;; Loop
  for jj=0L,nstrct-1 do begin
      ii = gd[jj]
      ;; For passing
      tmp = strct[ii]

      case strct[ii].name of 
          'Sun': begin
              tmp.param[0] = 5500. ;; Temperature
              outr_overlays_star, tmp, INIT=init, WVMNX=wvmnx, $
                                  FIDDLE=fiddle, PARAM=parm
          end
          'Star': outr_overlays_star, tmp, INIT=init, WVMNX=wvmnx, $
            FIDDLE=fiddle, PARAM=parm
          'Human': outr_overlays_star, tmp, INIT=init, PARAM=[300.], $
            WVMNX=wvmnx, FIDDLE=fiddle
          'Na Lamp': outr_overlays_lamp, tmp, 0
          'QSO': outr_overlays_qso,  tmp, INIT=init, $
            FIDDLE=fiddle, PARAM=parm
          else: stop
      endcase

      ;; Save
      strct[ii] = tmp
  endfor
  

  return
end
