;+ 
; NAME:
;  grb_fxlum
;   Version 1.1
;
; PURPOSE:
;    Calculates GRB flux as a function of time given the flux at
;    one point in time.  The solution is based on the paper by
;    Sari, Piran, & Narayan 1998.
;
; CALLING SEQUENCE:
;   
;   grb_fxlum, fini, tini, freq, RAD=, FLUX=, LUM=, Z=, $
;              tEarth=,  TEXP=, TSTART=, NSTP=
;
; INPUTS:
;     fini -- Flux measured at tini in muJy
;     tini -- Time at Earth when fini was measured (seconds)
;     freq -- Frequency of the radiation
; RETURNS:
;   flux= -- Array of flux values computed as function of time
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;    /RAD   -- Radiative expansion [Default: adiabatic]
;    z=     -- Redshift [default: 1.]
;  tstart=  -- Starting time of observation [default: 3600.]
;  tend=    -- Starting time of observation [default: 7200.]
;  nstp=    -- Number of time steps [default: 1000L]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   grb_fxlum, 1., 7200., 1e15
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   29-Oct-2003 Written by JXP (based on Sari, Piran, Narayan 1998)
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro grb_fxlum, fini, tini, freq, RAD=rad, FLUX=flux, LUM=lum, Z=z, $
               tEarth=tEarth,  TEXP=texp, TSTART=tstart, NSTP=nstp

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'grb_fxlum, fini, tini, freq, /RAD, FLUX=, z=, tEarth=, texp='
    print, '    TSTART=, NSTP= [v1.1]'
    return
  endif 

  ;; Optional keywords
  if not keyword_set( Z ) then z = 1.
  if not keyword_set( tstart ) then tstart = 3600.
  if not keyword_set( texp ) then texp = 7200.
  if not keyword_set( NSTP ) then nstp = 1000L

  flux = dblarr(nstp)
  tEarth = dblarr(nstp)

  ;; Adiabatic
  if not keyword_set( RAD ) then begin
      ;; Calculate timescales
      t_0 = 210.*24.*3600.  ; seconds
      t_m = 0.69 * (freq/1e15)^(-0.666666)*24.*3600.  ; seconds

      ;; Frequency 'cut'
      nu_0 = 1.8e11 
      if freq LT nu_0 then stop  ;; Too low for this calculation

      ;; Flux at t_m (muJy)
      if tini/(1.+z) LT t_m then fx_m = fini* (t_m*(1.+z)/tini)^(-0.25) $ 
      else fx_m = fini* (t_m*(1.+z)/tini)^(-1.375) 

      exp1 = 0.25
      exp2 = 1.375

  endif else begin  ;; Radiative

      ;; Calculate timescales
      t_0 = 4.6*24.*3600.  ; seconds
      t_m = 0.29 * (freq/1e15)^(-7./12.)*24.*3600.  ; seconds

      ;; Frequency 'cut'
      nu_0 = 8.5e12 
      if freq LT nu_0 then stop  ;; Too low for this calculation

      ;; Flux at t_m (muJy)
      if tini/(1.+z) LT t_m then fx_m = fini* (t_m*(1.+z)/tini)^(-0.571429) $ 
      else fx_m = fini* (t_m*(1.+z)/tini)^(-1.857) 

      exp1 = 4./7.
      exp2 = 1.857
  endelse

  ;; Loop on time to calculate flux
  for ii=0L,nstp-1 do begin
      ;; Time at GRB
      t = (tstart + texp*float(ii)/float(nstp)) / (1.+z)
      tEarth[ii] = t*(1.+z)
      
      ;; Flux
      if t LT t_m then flux[ii] = fx_m* (t_m/t)^exp1 $
      else flux[ii] = fx_m* (t_m/t)^exp2
      
  endfor


  ;; Luminosity
  dt = tEarth - shift(tEarth,1)
  dt[0] = dt[1]
  lum = total(flux*dt, /cumulative)

  return
end
