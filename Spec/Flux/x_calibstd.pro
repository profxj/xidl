;+ 
; NAME:
; x_calibstd   
;   Version 1.1
;
; PURPOSE:
;    Plots any array interactively
;
; CALLING SEQUENCE:
;   
;   spec = x_apall(ydat, [head])
;
; INPUTS:
;   ydat       - Values 
;   [head]     - Header
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   wave       - wavelength array
;   DISPLAY    - Display the sky subtracted image with xatv
;   OVR        - String array for ov region:  '[2050:2100, *]'
;   ERROR      - Variance array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_calibstd, kast, 0, 1, 0
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   29-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_calibstd, wave, flux, hstfil, outfil, CHKFIT=chkfit, EXP=exp, BSPLIN=bsplin,$
                GDWV=gdwv, BSET=bset, YFIT=yfit, SENS=sens

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'x_calibstd, wave, flux, hstfil, outfil, GDWV= [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( EVERYN ) then everyn = 5
  if not keyword_set( EXP ) then begin
      print, 'Assuming an exposure time of 1s'
      exp = 1.
  endif


; HST
  ;; Read file
  hst = xmrdfits(hstfil, 1, /silent)

  ;; BSpline
  if keyword_set( BSPLIN ) then begin
      bset = bspline_iterfit(hst.wavelength, hst.flux, yfit=yfit, everyn=everyn)
      if keyword_set( CHKFIT ) then begin
          x_splot, hst.wavelength, hst.flux, ytwo=yfit, /block
          stop
      endif
      ;; Calculate at wavelength
      sens = bspline_valu(wave, bset) / (flux / exp)
      
  endif else begin ;; SPLINE
      ;; Calculate at wavelength
      splin = spl_init(hst.wavelength, hst.flux, /double)
      std = spl_interp(hst.wavelength,hst.flux,splin, wave)
      if keyword_set( CHKFIT ) then $
        x_splot, hst.wavelength, hst.flux, xtwo=wave,ytwo=std, /block

      ;; Sens function
      sens = std / (flux / exp)
  endelse

  ;; Mask
  if keyword_set( GDWV ) then begin
      good = where( wave GT gdwv[0] AND wave Lt gdwv[1], ngood)
      if ngood EQ 0 then stop
      tmpwave = wave[good]
      sens = sens[good]
  endif else tmpwave = wave
  srt = sort(tmpwave)
  tmpwave = tmpwave[srt]
  sens = sens[srt]
  
  ;; BSpline
  bset = bspline_iterfit(tmpwave, sens, yfit=yfit, everyn=everyn)

  ;; Output
  mwrfits, bset, outfil, /create
  
  print, 'x_calibstd: All Done!'
  return
end
