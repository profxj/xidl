;+ 
; NAME: 
; igm_teff_limit   
;    Version 1.1
;
; PURPOSE:
;    Calculates tau_eff^LL for an emission redshift, an f(N) power-law
;    expression and a redshift (or set of redshifts)
;
; CALLING SEQUENCE:
;   
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /SUMZ_FIRST -- Sum on redshift first
;
; OPTIONAL OUTPUTS:
;  teff_limit -- Array of tau_eff^LL values from z912 to zem (or in
;                NHI)
;
; COMMENTS:
;
; EXAMPLES:
;   
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   June-2011 Written by JXP
;-
;;  Calculates the effecitive opacity from Lyman limit absorption
;;;;;
;;;;
function igm_teff_limit, z912, zem, LYMAN=lyman, POWERFN=powerfn, N_eval=n_eval, $
                         TEFF_LIMIT=teff_limit, FNZ=fnz, DEBUG=debug, ZVAL=zval, $
                         SUMZ_FIRST=sumz_first, NVAL=lgNval, W05MAP=w05map, verbose=verbose, $
                         H0=h0
  
  if not keyword_set(N_eval) then N_eval = 1000L
  c = x_constants()

  if not keyword_set(POWERFN) then stop; init_madau95_fnz, powerfn

  
  ;; NHI array
  lgNval = 11.5 + 10.5*dindgen(N_eval)/(N_eval-1) ;; This is base 10 [Max at 22]
  dlgN = lgNval[1]-lgNval[0]
  Nval = 10.d^lgNval

  ;; z array
  zval = z912 + (zem-z912)*findgen(N_eval)/(N_eval-1)
  dz = abs(zval[1]-zval[0])
  
  teff_limit = dblarr(N_eval)

  dXdz = cosm_dxdz(zval, silent=(not keyword_set(VERBOSE)), W05map=w05map, H0=H0)
  if keyword_set(FNZ) then dXdz = replicate(1.,N_eval)

  ;; Evaluate f(N,X)
  velo = (zval-zem)/(1+zem) * (c.c/1e5) ; Kludge for eval_powerfn [km/s]
  log_fnX = eval_powerfn(POWERFN, lgNval, zem, VEL_ARRAY=velo)  
  log_fnz = log_fnX + replicate(1., N_eval) # alog10(dxdz)
  
  ;; Evaluate tau(z,N)
  teff_engy = c.Ryd/c.eV
  teff_cross = x_photocross(1,1,teff_engy)
  sigma_z = teff_cross * ((1+zval)/(1+zem))^(2.75)  ;; Not exact but close
  tau_zN = Nval # sigma_z  

  ;; Integrand
  intg = 10.d^(log_fnz) * (1. - exp(-1.*tau_zN))

  if not keyword_set(SUMZ_FIRST) then begin
     ;; Sum in N first
     N_summed = total(intg*(Nval # replicate(1.,N_eval)),  1) * dlgN * alog(10.)
     
     ;; Sum in z
     teff_limit = reverse(total(reverse(N_summed),/cumul)) * dz 
  endif else begin
     z_summed = total(intg, 2) * dz
     teff_limit = total(reverse(z_summed * Nval), /cumul) * dlgN * alog(10.)
  endelse
  ;stop

  ;; Debug
  if keyword_set(DEBUG) then begin
;        x_splot, lgNval, alog10(10.d^(log_fnX) * dxdz * dz * Nval), /bloc
;        x_splot, lgNval, total(10.d^(log_fnX) * dxdz * dz * Nval,/cumul) * dlgN * alog(10.) / teff_lyman[qq], /bloc
     printcol, lgnval, log_fnx, dz,  alog10(10.d^(log_fnX) * dxdz * dz * Nval)
     writecol, 'debug_file'+strtrim(qq,2)+'.dat', $
               lgNval, restEW, log_fnX
     stop
  endif
     
  return, max(teff_limit)
end
