;+ 
; NAME:
; x_fluxcalib   
;    Version 1.1
;
; PURPOSE:
;    Given a spectrum and a flux calib solution, calibrate
;
; CALLING SEQUENCE:
;   fx_fnu = x_fluxcalib(wv, fx, fitstr, [var, newvar], /FLAMBDA, DLMB=, TRUCONV=)
;
; INPUTS:
;   wv   - Wavelength array
;   fx   - Stellar flux  (e- per pixel)
;   fitstr - Calibration fitting function (assumes alog10 for flux)
;   [var]  - Variance array
;
; RETURNS:
;   fx_fnu - Flux in fnu (or flambda)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /FLAMBDA - Return flambda instead of fnu
;   DLMB -- Delta lambda (or km/s if negative) of wavelength array
;   /LINEAR - Assume fitstr *not* in alog10 for flux
;   
; OPTIONAL OUTPUTS:
;   newsig -- Fluxed variance array
;
; COMMENTS:
;
; EXAMPLES:
;   fx_fnu = x_fluxcalib(wv, fx, fitstr)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_fluxcalib, wv, fx, fitstr, var, newsig, LINEAR=linear, $
                      DLMB=dlmb, FLAMBDA=flambda, TRUCONV=truconv
;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'fx_fnu = x_fluxcalib(wv, fx, fitstr, [var, newvar], /FLAMBDA, '
    print, '            /LINEAR, DLMB=, TRUCONV=) [v1.1]'
    return, -1
  endif 

;  Optional Keywords
  if not keyword_set(scale) then scale = 1.

; Calculate the conversion factor
  convfact = x_calcfit(wv, FITSTR=fitstr)
  ; Take 10^
  if keyword_set(LINEAR) then truconv = convfact $
  else truconv = 10.^convfact   ; could be really large

; Find dwv at each pixel

  npix = n_elements(wv)
  if keyword_set(DLMB) then begin
      if dlmb GT 0 then dwv = replicate(dlmb, npix) $  ;  Constant dlamb
        else dwv = -(dlmb/2.9979e5)*wv         ; Constant delv (km/s)
  endif else begin
      dwv = (shift(wv,-1)-shift(wv, 1))/2.
      dwv[0] = wv[1]-wv[0]
      dwv[npix-1] = wv[npix-1]-wv[npix-2]
  endelse
      
; Convert obj to fnu
  obj_fnu = (fx/dwv) * wv^2 / (3e18)
  if not keyword_set( FLAMBDA ) then fin_flux = obj_fnu*truconv $  ; fnu
  else fin_flux = obj_fnu*truconv*3.e18/(wv^2)   ; flambda

; VARIANCE
  if keyword_set( VAR ) then begin
      sig = var
      a = where(var GT 0.)
      sig[a] = sqrt(var[a])
      obj_sig = (sig/dwv) * wv^2 / (3e18)
      if not keyword_set( FLAMBDA ) then newsig = obj_sig*truconv $ ; fnu
      else newsig = obj_sig*truconv*3.e18/(wv^2) ; flambda
      badsig = where(newsig LT 0., nbad)
      if nbad NE 0 then newsig[badsig] = -1.
  endif

; Return
  return, fin_flux


      
end
      
