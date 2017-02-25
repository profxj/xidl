;+
;
;
; For a given halo mass, redshift and  cosmology, find an estimate of virial radius,
; velocity following Bryan and Norman 1998
;
; halo_mass         halo mass in Msun
; red               redshift 
; EXTRA             cosmology
; H0                hubble constant
;
;
;-


pro m_cosm_virialrel, halo_mass, red, _EXTRA=EXTRA, H0=H0, silent=silent

  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, cosm_L, cosm_r, sigma_8
  
  ;;set intialize cosmology
  cosm_common, H0=h0, _EXTRA=extra, SILENT=silent

  
  ;;derive expansion factor
  aexp=1./(1+red)
  splog, 'Expansion factor ', aexp

  ;;find cosmology relation 
  little_h=(cosm_h/100.)
  Ez_square=cosm_dm*(1+red)^3+cosm_L ;;neglect Omega_r
  Hz=cosm_h*sqrt(Ez_square)                ;;km/s/Mpc
  Omega_mz=cosm_dm*(1+red)^3/Ez_square
  x=Omega_mz-1
  Delta_crit=18.*!pi^2+82*x-39*x^2
  
  splog, 'Virial overdensity ', Delta_crit 
  splog, 'H(z) ', Hz
  splog, 'Omega_m(z) ', Omega_mz
  
  ;;find critical density
  rho_crit=1.88D-29*little_h^2*Ez_square   ;;g/cm^3
  splog, 'Critical mass density (g/cm^3) ', rho_crit


  ;;find virial radius
  halo_mass_g=halo_mass*1.98892D33 ;;gram
  rvir_cm=((3*halo_mass_g)/(4*!pi*Delta_crit*rho_crit))^(1./3)  ;;cm 
  rvir=rvir_cm*3.24077649D-22  ;;kpc
  splog, 'Virial radius (kpc) ', rvir

  ;;find virial velocity
  v_circ=sqrt(6.6726D-8*halo_mass_g/rvir_cm)*1D-5 ;;km/s
  splog, 'Virial velocity (km/s) ', v_circ


  
end

