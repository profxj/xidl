;+ 
; NAME:
; ism_pltlin
;   Version 1.0
;
; PURPOSE:
;    Plots an ISM transition (e.g. SiII 1526)
;
; CALLING SEQUENCE:
;  ism_pltlin, lambda, N, b, XMNX=, YMNX=, EW=, /NOPLOT
;
; INPUTS:
;   lambda -- Rest wavelength of the transition
;   N -- Column density [log]
;   b -- Doppler parameter [km/s]
;
; RETURNS:
;  Plots to screen
;   
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  XMNX=  -- Set plot x range
;  YMNX=  -- Set plot y range
;  VMNX=  -- Set plot x range using velocities [km/s]
;
; OPTIONAL OUTPUTS:
;  EW=  -- EW of the transition
;
; COMMENTS:
;
; EXAMPLES:
;   cldy_Uplot, grid, 19.0d, -1.0d, -1.0d, [[14,2], [14,3]]
;
;
; PROCEDURES/FUNCTIONS CALLED:
; getcolor
; x_setclrs
; getabnd
;
; REVISION HISTORY:
;   15-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ism_pltlin, lambda, N, b, XMNX=xmnx, YMNX=ymnx, $
                EW=ew, NOPLOT=noplot, VMNX=vmnx, NPIX=npix, $
                DLAMB=dlamb

;
  if  N_params() LT 3  then begin 
      print, 'Syntax - ' +$
        'ism_pltlin, lambda, N, b, XMNX=, YMNX=, EW=, VMNX=, /NOPLOT [v1.1]'
      return
  endif 

; Optional keywords

  if not keyword_set( YMNX ) then ymnx = [-0.1, 1.05]
  if not keyword_set(NPIX) then npix = 5000
  if not keyword_set(DLAMB) then dlamb = 10

  lin = x_setline(lambda, /close)
  lin.N = N
  lin.b = b

  print, 'Atomic: ', lin.wrest, lin.f

  ;; Create wavelength array
  wv = (lambda - dlamb) * 10^(dindgen(npix)*1.883e-6)
  fx = x_voigt(wv, lin, FWHM=4.)

  if keyword_set(VMNX) then begin
      vel = (wv-lambda)/lambda * 3e5
      mn = min(abs(vel-min(vmnx,max=mxv)),i1)
      mn = min(abs(vel-mxv),i2)
      xmnx = [wv[i1], wv[i2]]
  endif


  if not keyword_set(XMNX) then begin
      dwv = lin.wrest* b / 3e5
      xmnx = lin.wrest + 5*dwv*[-1,1]
  endif

  clr = getcolor(/load)

  if not keyword_set(NOPLOT) then begin
      plot, wv, fx, psym=10, xrange=xmnx, thick=3, $
            yrange=ymnx, color=clr.black, charsiz=1.3, $
            background=clr.white, ystyle=1, xstyle=1, $
            xtitle='Wavelength', ytitle='Normalized Flux'

      xyouts, xmnx[0]+0.1*(xmnx[1]-xmnx[0]), 0.1, lin.ion, color=clr.black, $
              charsiz=1.5
      oplot, xmnx, [0., 0.], color=clr.green, linest=2, thick=4
  endif

  ;; EW
  dwv = wv - shift(wv,1)
  dwv[0] = dwv[1]
  ew = total( (1.-fx)*dwv)
  print, 'ism_pltlin:  EW = ', ew, ' Ang'

  return
end

