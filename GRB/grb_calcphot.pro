;+ 
; NAME:
;  grb_calcphot
;   Version 1.1
;
; PURPOSE:
;  Calculate the number of photons emitted in an energy interval dE
;  during a time interval dt.  The user inputs a structure which
;  describes the GRB afterglow light curve.
;
; CALLING SEQUENCE:
;   
;   total_phot = grb_calcphot( dE, dt, grb_after, Rphot=, nphot=)
;
; INPUTS:
;     dE  -- Energy interval (GRB frame)  [eV]   
;     dt  -- time interval (observer frame)
;     grb_after -- Afterglow structure
;
; RETURNS:
;   total_phot -- Number of photons emitted in the energy and time
;                 intervals
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  Rphot= -- Radius at which to calculate the column density of photons
;  nphot= -- Column density of photons at pecified distance
;  /SILENT -- Suppress cosmo info
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Feb-2006 Written by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function grb_calcphot, dE, dt, grb_after, Rphot=Rphot, nphot=nphot, SILENT=silent

common grb_splinelc, grbsp_nu, grbsp_mag, grbsp_t, grbsp_splin

  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'nphot = grb_calcphot( dE, dt, grb, Rphot=, /SILENT) [v1.0]'
    return, -1
  endif 

;  cd, getenv('XIDL_DIR')+'/GRB/Light_Curves/', CURRENT=current
;  RESOLVE_ROUTINE, 'grb_splnlc', /COMPILE_FULL_FILE, /EITHER
;  cd, current
;  RESOLVE_ALL, RESOLVE_FUNCTION='grb_splnlc_spline'


;  dE = [7.7, 9.54]
;  dt = [50., 3840]
;  grb = {grbafterglow}
;  grb.beta = -0.5
;  grb.t0 = 50
;  grb.alpha = -0.8
;  grb.mag_wv = 6588.
;  grb.mag = 13.5  ; AB
;  grb.z = 1.5495
;  grb_after = grb

  ;; Optional keywords
  if not keyword_set(h0) then h0 = 70.
  
  if n_elements(dt) NE 2 then stop  ;; Had some funny overloading going on
  c = x_constants()
      
  ;; Convert Energy to frequency
  dnu = (dE*c.eV)/c.h / (grb_after.z + 1) ; Ev -> Hz and Redshift

  nu_mag = c.c / (grb_after.mag_wv * 1e-8)  ; Hz
  

  ;; Power law
  ;; Fnu = F0 (nu_0/nu)^-beta (t/t_0)^alpha
  if grb_after.alpha LT 9.99 then begin
      ;; Convert AB magnitude to Flux at Earth
      flux = 10^( -1.*(grb_after.mag + 48.6)/2.5)

      ;; Integrate 
      const = flux / (grb_after.t0^(-1.*grb_after.alpha)) / c.h / $
              (nu_mag^(-1.*grb_after.beta)) / $
              (-1.*grb_after.beta) / (1-grb_after.alpha)
      const = abs(const)
      
      ;; Assumer power law
      integral = abs( dnu[0]^(-1.*grb_after.beta) - dnu[1]^(-1.*grb_after.beta)) $
                 * abs( dt[0]^(1-grb_after.alpha) - $
                        dt[1]^(1-grb_after.alpha))
  endif else begin ;; Spline for time, but (nu_0/nu)^-beta for frequency

      ;; Integrate 
      const = 1. / c.h / (nu_mag^(-1.*grb_after.beta)) / $
              (-1.*grb_after.beta) 
      const = abs(const)

      ;; Time interval (log steps)
      nstp = 10000L
      tarr = dt[0] * $
             exp(alog(dt[1]/dt[0]) * dindgen(nstp) / nstp )

      ;; Get magnitudes
      if not keyword_set(grbsp_t) then grb_splnlc, grb_after
      marr = grb_splnlc_spline(tarr)

      ;; Convert to flux and add it up
      flux = 10^( -1.*(marr + 48.6)/2.5)
      tot_flux = int_tabulated(tarr, flux, /double)
      
      ;; Assumer power law
      integral = abs( dnu[0]^(-1.*grb_after.beta) - $
                      dnu[1]^(-1.*grb_after.beta)) $
                 * tot_flux
  endelse

  nphot = const * integral  ; cm^-2

  lum_dist = cosm_dist(grb_after.z, H0=h0, /init, /Lum, SILENT=silent)  ; Mpc
  tot_phot = 4.d * !pi * (lum_dist*c.Mpc)^2 * const * integral $
             / (1.+grb_after.z)^2 ; # photons

  ;; Rphot
  if keyword_set(Rphot) then begin
      if Rphot GT 0. then nphot = tot_phot / (4 * !pi * (Rphot * c.pc)^2) ; cm^-2
      print, 'Photons at ', Rphot, 'pc are ', nphot, ' cm^-2' 
  endif
      
  return, tot_phot
end
