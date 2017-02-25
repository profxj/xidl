;+ 
; NAME:
; dla_si2star   
;   Version 1.1
;
; PURPOSE:
;    Given a normalized QSO spectrum (fits file) and the redshift of
;    the DLA, overplots a SiII* feature at the expected spot based on
;    the observed CII* profile.  THe program also plots the 
;    CII* profile
;
; CALLING SEQUENCE:
;  dla_si2star, fits_fil, zabs, VMNX=, PSFIL=, YMNX=, RTIO=
;
; INPUTS:
;  fits_fil -- FITS file containing the QSO data [Assumes HIRES
;              format]
;  zabs  -- Absorption redshift
;
; RETURNS:
;
; OUTPUTS:
;   PSFIL=  -- Writes PS file 
;
; OPTIONAL KEYWORDS:
;  RTIO= -- Ratio of optical depth of SiII* to CII* [deafult: 0.05]
;  VMNX= -- Velcoity region to plot [default: -100 to 100 km/s]
;  YMNX= -- Ymin, ymax of plot [deafult: 0., 1.]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   dla_si2star, 'Q1331_f.fits', 1.770
;
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   29-May-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dla_si2star, fits_fil, zabs, VMNX=vmnx, PSFIL=psfil, YMNX=ymnx, RTIO=rtio

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'dla_si2star, fits_fil, zabs, VMNX=, PSFIL=, YMNX=, RTIO=  [v1.1]'
    return
  endif 

; Optional Keywords
  if not keyword_set( VMNX ) then vmnx = [-100.d, 100.]
  if not keyword_set( YMNX ) then ymnx = [[0.0, 1.1], [0.0,1.1]]
  if not keyword_set( RTIO ) then rtio = 0.05

  error_fil = strmid(fits_fil, 0L, strlen(fits_fil)-6)+'e.fits'
; Read in data
  fx = x_readspec(fits_fil, FIL_SIG=error_fil, SIG=sig, WAV=wave)

  clr = getcolor(/load)
  !P.MULTI= [0,1,2,0,1]

  if keyword_set(PSFIL) then begin
      device, get_decomposed=svdecomp
      device, decompose=0
      ps_open, file=psfil, font=1, /color
  endif

; Create velocity arrays
  wrest = [1264.7377d, 1335.7077d]

  vel_arr = x_allvelo(wave, zabs, wrest, vmnx, all_pmnx=all_pmnx)

  spaces = replicate('!17 ',30)
; Plot
  plot, vel_arr[0:all_pmnx[2,0],0], $
    fx[all_pmnx[0,0]:all_pmnx[1,0]], xrange=vmnx, $
    yrange=ymnx[*,0], xtickn=spaces, xmargin=[9,3], $
    ymargin=[0,0], $
    charsize=1.8, psym=10, background=clr.white, color=clr.black, $
    xstyle=1, ystyle=1

  oplot, [-999.,999.], [1.,1.], color=clr.red, linestyle=2
  xyouts, vmnx[0]+5., ymnx[0,0]+0.1, 'SiII*', charsize=3.

; Overplot fake SiII*
  a = where(fx[all_pmnx[0,1]:all_pmnx[1,1]] LT 0.8, na)
  if na NE 0 then begin
      oplot, vel_arr[a], exp(rtio*alog(fx[all_pmnx[0,1]+a]>0.03)), $
        color=clr.blue, psym=-1
  endif

  plot, vel_arr[0:all_pmnx[2,1],0], $
    fx[all_pmnx[0,1]:all_pmnx[1,1]], xrange=vmnx, $
    yrange=ymnx[*,1], xmargin=[9,3],  ymargin=[3,0], $
    charsize=1.8, psym=10, background=clr.white, color=clr.black, $
    xstyle=1, ystyle=1

  oplot, [-999.,999.], [1.,1.], color=clr.red, linestyle=2
  xyouts, vmnx[0]+5., ymnx[0,1]+0.1, 'CII*', charsize=3.

  !P.MULTI= [0,1,1,0,1]
  if keyword_set( PSFIL ) then begin
      ps_close, /noprint, /noid
      device, decomposed=svdecomp
  endif

; Calculate optical depths

  ;; CII*
  tau_CII = -alog(fx[all_pmnx[0,1]:all_pmnx[1,1]])
  tot_tCII = total(tau_CII)

  var_tau = (sig[all_pmnx[0,1]:all_pmnx[1,1]] / fx[all_pmnx[0,1]:all_pmnx[1,1]])^2
  tot_stCII = sqrt(total(var_tau))

  print, 'tau(CII*): ', tot_tCII, tot_stCII

  ;; SiII*
  gd = where(sig[all_pmnx[0,0]:all_pmnx[1,0]] GT 0.)
  tau_SiII = -alog(fx[all_pmnx[0,0]+gd])
  tot_tSiII = total(tau_SiII)

  var_tau = (sig[all_pmnx[0,0]+gd] / fx[all_pmnx[0,0]+gd])^2
  tot_stSiII = sqrt(total(var_tau))
  print, 'tau(SiII*): ', tot_tSiII, tot_stSiII

  ;; Nsig
  print, 'Nsig = ', (tot_tCII*rtio-(tot_tSiII>0.))/tot_stSiII

  return
end
