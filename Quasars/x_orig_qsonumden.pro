;+ 
; NAME:
; x_qsonumden
;
; PURPOSE:
;    Gives #qsos per comoving Mpc-3 at a given redshift for 
;    a given magnitude limit (at 1450A).  Assumes a cosmology
;    H0 = 65 km/s/Mpc; Omega = 0.35, Lambda = 0.65 (following Fan et
;    al. 2004).
;    psi(L) = psistar * L^(beta)   :: beta = -3.2 [default]
;    The luminosity function is assumed to be a double power-law
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;
; RETURNS:
;  NUM_DEN -- Number density of quasars per Mpc^3 at z to Mlim
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  ZEVOL=   ;; Exponential evolution of L* [default: -0.43]
;  PSI267 = ;; Number density of QSOs at z=6 for M<=-26.7
;              [default:6d-10]
;  BETA=    ;; Slope of the bright end [default: -3.2]
;  MSTAR=   ;; Value for break in power-law (At 1450A) 
;              [default: -24.5]
;  BLOW=    ;; Power-law exponenent for L < L_*
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   13-Feb-2008 Written by JXP
;-
;------------------------------------------------------------------------------
function x_qsonumden, z, Mlim, PSISTAR=psistar, BETA=beta, ZEVOL=zevol, $
  PSI267=PSI267, MSTAR=Mstar, BLOW=blow, NOCOSM=nocosm

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'numqso = x_qsonumden(z, Mlim) [v1.0]'
    return, -1
  endif 

  if not keyword_set(BLOW) then blow = -1.64  ;; Power-law at low luminosity
  if not keyword_set(BETA) then beta = -3.2
  if not keyword_set(ZEVOL) then zevol = -0.43
  if not keyword_set(MSTAR) then MSTAR = -24.5 ;; (at 1450A)
  if not keyword_set(PSI267) then psi267 = 6d-10
     ;; Mpc^-3, H0 = 65 km/s/Mpc; Omega = 0.35, Lambda = 0.65
     ;; Fan et al. 2004

  c = x_constants()
  num_den = dblarr(n_elements(Mlim))

  stop ;; JXP -- This code is likely buggy.  Best to avoid  (Nov 2011)

  ;; Luminosity distance
  if not keyword_set(NOCOSM) then $
    cosm_common, H0=65., OMegaDM = 0.35, Omegavac=0.65, /silent
  r1 = cosm_dist(z, /silent)
  dL = r1 * (1+z)               ; Luminosity distance Mpc

  ;; Find psistar
  psistar = -1.*(beta+1) * PSI267 * (10.^(-0.4*(Mstar + 26.7)))^(beta+1)

  ;; Find apparent MSTAR
  dum = x_abtoflam(MSTAR, 1450., FNU=f1450)
  L1450 = 4*!dpi*(10*c.pc)^2 * f1450
  f_filter = L1450 / (4*!dpi*(dL*c.mpc)^2) * (1+z)
  AMSTAR = -2.5 * alog10(f_filter) - 48.6  ;; At 1450

  ;; Integrate below M*
  low = where(Mlim GT AMSTAR, nlow)
  if nlow GT 0 then begin
      num_den[low] = -1. * psistar * $
                     (1. - (10.^(0.4*(Mlim[low] - AMSTAR)))^(blow+1)) $
                     / (blow+1)
  endif
;  print, num_den

  ;; Set integral limit for bright end
  MHIGH = Mlim < AMSTAR

  ;; Integrate above M*
  num_den = num_den - 1. * psistar * (10.^(-0.4*(MHIGH - AMSTAR)))^(beta+1) $
            / (beta+1)
;  stop
;  print, num_den

  
  lum_evol = exp(ZEVOL* (z-6)) ;; Modulate PSISTAR by redshift 
  
;; Number density
  
  num_den = num_den * LUM_EVOL

  return, num_den
end
