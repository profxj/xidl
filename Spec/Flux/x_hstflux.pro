;+ 
; NAME:
; x_hstflux   
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

function x_hstflux, wv, fx, fitstr, var, newsig, $
                      DLMB=dlmb, FNU=fnu, TRUCONV=truconv, EXP=exp
;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'fx_fnu = x_fluxcalib(wv, fx, fitstr, [var, newvar], /FNU, '
    print, '            DLMB=, TRUCONV=) [v1.1]'
    return, -1
  endif 

;  Optional Keywords
  if not keyword_set( EXP ) then exp = 1.

;; THIS PROGRAM USED TO ASSUME THE DWV OF THE INPUT IS THE SAME AS THE FIT
  
  ;; Dwv  (Ang)
  npix = n_elements(wv)
  if wv[1] GT wv[0] then begin
      dwv = (shift(wv,-1)-shift(wv, 1))/2.
      dwv[0] = wv[1]-wv[0]
      dwv[npix-1] = wv[npix-1]-wv[npix-2]
  endif else begin
      dwv = (shift(wv,1)-shift(wv, -1))/2.
      dwv[0] = wv[0]-wv[1]
      dwv[npix-1] = wv[npix-2]-wv[npix-1]
  endelse

; Calculate the conversion factor
  convfact = bspline_valu(wv, fitstr)

; Convert obj to fnu
  obj_flmb = fx * convfact / dwv / exp
  if keyword_set( FNU ) then fin_flux = obj_flmb/3.e18*(wv^2) $
  else fin_flux = obj_flmb     ; flambda

  ;; VARIANCE
  if keyword_set( VAR ) then begin
      sig = var
      a = where(var GT 0.)
      sig[a] = sqrt(var[a])
      obj_sig = sig * convfact / exp
      if not keyword_set( FNU ) then newsig = obj_sig $
      else newsig = obj_sig/3.e18*(wv^2) ; fnu
      badsig = where(newsig LT 0., nbad)
      if nbad NE 0 then newsig[badsig] = -1.
  endif

  ;; Return
  return, fin_flux


      
end
      
