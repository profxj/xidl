;+ 
; NAME:
; cosm_common
;
; PURPOSE:
;    Routine to initialize and set values in the Cosmology common
;    block named 'cosmolgy_cmmn'
;
; CALLING SEQUENCE:
;   cosm_common
;
; INPUTS:
;   H0 =   Hubbles constant in km/s/Mpc
;   Omegavac = Lambda value   [Default: 0.7]
;   OmegaDM = Omega for Dark Matter  [Default: 0.3]
;
; RETURNS:
;
; OUTPUTS:
;   
; OPTIONAL KEYWORDS:
;  /W05MAP -- Applies the 5yr WMAP cosmolgy (Dunkeley et al. 2009) 
;  /SILENT -- No screen output
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
;   22-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro cosm_common, H0=h0, Omegavac=omegavac, OmegaDM=omegaDM, $
                 W05MAP=W05MAP,  SILENT=silent, VANILLA=vanilla

  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, $
     cosm_L, cosm_r, sigma_8

;; Original defaults are 0.3, 0.7 and 75km/s/Mpc
  if not keyword_set( OmegaDM ) then cosm_dm = 0.3 else cosm_dm = omegaDM
  if size( H0 , /type) EQ 0 then cosm_h = 75. else cosm_h = H0
  if size( Omegavac, /type) EQ 0 then cosm_L = 0.7 else cosm_L = Omegavac

  if keyword_set(W05MAP) then begin  ;; 5yr WMAP;  Dunkley et al. 2009
      cosm_h = 72. 
      if keyword_set(H0) then cosm_h = H0  ;; For QPQ6 in particular
      cosm_dm = 0.26
      cosm_L = 0.74
      sigma_8 = 0.796
      omega_Ob = 0.0455
   endif

  if keyword_set(VANILLA) then begin  
      cosm_h = 70. 
      cosm_dm = 0.30
      cosm_L = 0.70
      sigma_8 = 0.80
      omega_Ob = 0.0455
   endif

  cosm_K = 1.d - cosm_dm - cosm_L
  cosm_r = 2.4725753d-05

  if not keyword_set(SILENT) then begin
      print, 'cosm_common: Using this cosmology --'
      print, 'Omega_m = ', cosm_dm
      print, 'H (km/s/Mpc) = ', cosm_h
      print, 'Omega_L  = ', cosm_L
  endif
  return

end

