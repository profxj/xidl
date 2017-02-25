;+ 
; NAME:
;  grb_calclum
;   Version 1.1
;
; PURPOSE:
;    Given a grb structure (magnitude at a time t,  z, spectral slope)
;    calculates the luminsoity at any time  (erg/s/cm^2/Hz)
;
; CALLING SEQUENCE:
;   
;   lum_nu = grb_calclum(grb, nu, [t])
;
; INPUTS:
;     grb -- GRB structure
;     nu  -- Frequency (GRB frame)
;     [t] -- Time (observer frame). Default = t0 in GRB structure
;
; RETURNS:
;   lum_nu= 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; H0 -- Set Hubble constant
; TAVG=  -- Time interval to calculate average luminosity across
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

function grb_calclum, grb, nu, t, H0=h0, TAVG=tavg

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'lum_nu = grb_calclum(grb, nu, [t_obs], H0=, TAVG=) [v1.0]'
    return, -1
  endif 

  ;; Optional keywords
  if not keyword_set(t) then t = grb.t0
  if not keyword_set(h0) then h0 = 70.
  c = x_constants()
      
  ;; Convert magnitude to flux (assume AB)
  nu_mag = c.c / (grb.mag_wv * 1d-8)  ; Hz
  nuE = nu / (1. + grb.z)

  fnu = 10^( -1.*(grb.mag + 48.6)/2.5) * $
         (nuE/nu_mag)^(-1.*grb.beta) * $
         (t[0]/grb.t0)^(-1.*grb.alpha) ; ergs/s/cm^2/Hz

  ;; Time average?
  if keyword_set(TAVG) then begin
      if n_elements(t) NE 2 then stop
      fact = (t[1]^(1-grb.alpha) - t[0]^(1-grb.alpha)) / $
             (1-grb.alpha) / (t[1]-t[0]) / grb.t0^(-1.*grb.alpha) ;; Average
      fact = abs(fact) / $
             (t[0]/grb.t0)^(-1.*grb.alpha)  ;; Remove time factor from above
      ;; Result
      fnu = fnu * fact
  endif
  
  ;; Convert to Luminosity
  lum_dist = cosm_dist(grb.z, H0=h0, /init, /Lum)  ; Mpc

  ;; 1s at Earth is 1s in the GRB frame
  lum_nu = fnu * 4 * !pi * (lum_dist*c.mpc)^2 / (grb.z + 1.)  

  return, lum_nu
end
