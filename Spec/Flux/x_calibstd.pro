;+ 
; NAME:
; x_calibstd   
;   Version 1.1
;
; PURPOSE:
;  Given a standard star spectrum and a calibration file, calculate a
;  sensitivity function parameterized by a BSLINE or some other function.
;
; CALLING SEQUENCE:
;  x_calibstd, wave, flux, outfil, HSTFIL=, /CHKFIT, 
;               EXP=, BSPLIN=, SWV=, SFX=, 
;               GDWV=, BSET=, YFIT=, SENS=, EVERYN=
;
; INPUTS:
;  wave -- Wavelength array of standard star
;  flux -- Flux array of standard star
;
; RETURNS:
;
; OUTPUTS:
;  outfil -- FITS file to write sensitivity function 
;
; OPTIONAL KEYWORDS:
;  EXP= -- Exposure time [default: 1.]
;  HSTFIL= -- HST calibration file
;  /CHKFIT -- Plot the fit to the sensitivity function
;  EVERYN= -- Spacing of b-spline breakpoints in pixel space [default: 5]
;
; OPTIONAL OUTPUTS:
;  SWV=  -- Wavelength array of calibration data
;  SFX=  -- Flux array of calibration data
;  SENS= -- Sensitivity function
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

pro x_calibstd, wave, flux, outfil, HSTFIL=hstfil, CHKFIT=chkfit, $
                EXP=exp, BSPLIN=bsplin, SWV=swv, SFX=sfx, $
                GDWV=gdwv, BSET=bset, YFIT=yfit, SENS=sens, EVERYN=everyn,$
                ESO_FIL=eso_fil, MAB=mab, INTERFIT=interfit, NORD=nord

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'x_calibstd, wave, flux, outfil, hstfil=, /MAB, /CHKFIT, GDWV=, ESO_FIL= [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( EVERYN ) then everyn = 50
  if not keyword_set( EXP ) then begin
      print, 'Assuming an exposure time of 1s'
      exp = 1.
  endif

  ;; Dwv  (Ang)
  npix = n_elements(wave)
  if wave[1] GT wave[0] then begin
      dwv = (shift(wave,-1)-shift(wave, 1))/2.
      dwv[0] = wave[1]-wave[0]
      dwv[npix-1] = wave[npix-1]-wave[npix-2]
  endif else begin
      dwv = (shift(wave,1)-shift(wave, -1))/2.
      dwv[0] = wave[0]-wave[1]
      dwv[npix-1] = wave[npix-2]-wave[npix-1]
  endelse

; HST
  ;; Read file
  if keyword_set( HSTFIL ) then begin
      hst = xmrdfits(hstfil, 1, /silent)
      swv = hst.wavelength
      sfx = hst.flux
  endif

  ;; ESO/HST format
  if keyword_set(ESO_FIL) then $
    readcol, eso_fil, swv, sfx, format='D,F'

  ;; Mab
  if keyword_set(MAB) then begin
      fnu = 10^(-1.*(sfx + 48.6)/2.5)
      c = x_constants()
      sfx = fnu * c.c / swv^2 * 1e8
  endif
  

  ;; BSpline
  if keyword_set( BSPLIN ) then begin
      bset = bspline_iterfit(swv, sfx, yfit=yfit, everyn=everyn)
      if keyword_set( CHKFIT ) then begin
          x_splot, swv, sfx, ytwo=yfit, /block
          stop
      endif
      ;; Calculate at wavelength
      sens = bspline_valu(wave, bset) / (flux / exp / dwv)
      stop
      
  endif else begin ;; SPLINE
      ;; Calculate at wavelength
      splin = spl_init(swv, sfx, /double)
      std = spl_interp(swv, sfx, splin, wave)

      ;; Sens function
      sens = std / (flux / exp / dwv)  
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
  ivar = fltarr(n_elements(sens))
  ivar[*] = 1.

  ;; Abs
  bd = where((tmpwave GT 6864. AND tmpwave LT 6920.) OR $
             (tmpwave GT 7585. AND tmpwave LT 7684.), nbd)
  if nbd NE 0 then ivar[bd] = -1.

  if not keyword_set(INTERFIT) then begin
      bset = bspline_iterfit(tmpwave, sens, yfit=yfit, everyn=everyn, $
                             lower=2.5, upper=2.5, invvar=ivar, nbpkts=nord)
  endif else begin
      bset = bspline_iterfit(tmpwave, sens, yfit=yfit, nbkpts=nord)
;      yfit = x1dfit(tmpwave, sens, /inter, fitstr=fitstr)
;      bset = *fitstr.ffit
  endelse

  if keyword_set( CHKFIT ) then begin
      x_splot, tmpwave, sens, ytwo=yfit, /block
      stop
  endif

  ;; Output
  mwrfits, bset, outfil, /create
  
  print, 'x_calibstd: All Done!'
  return
end
