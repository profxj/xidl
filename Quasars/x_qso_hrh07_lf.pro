;+ 
; NAME:
; x_qso_hrh07_lf
;
; PURPOSE:
;    Returns the luminosity function of quasars for a default (or
;    specified) range of luminosity or magnitude at a given redshift.
;    Uses the "Full" model of Hopkins, Richards & Hernquist 2007. 
;    Tested against Phil's qlf_calculator.c code.
;
; CALLING SEQUENCE:
;  phi = x_qso_hrh07_lf(z, LVAL=)
;
; INPUTS:
;  z -- Redshift for evaluation
;
; RETURNS:
;  phi(L or M)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
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
;   Nov-2011 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_qso_hrh07_lf_lband, log_l_bol, mode, NU_EFF=nu_eff

  x = log_l_bol - 10.         
  lband = 0.                    ;

  case MODE of 
     0: begin
        lband = 1.
        nu_eff = 1.1992000e15 
     end
     -1: begin
        P0=8.99833 & P1=6.24800 & P2=-0.370587 & P3=-0.0115970;}
        nu_eff = 6.8136364e14 ;} // 4400 angstrom
     end
     -2: begin
        P0=10.6615 & P1=7.40282 & P2=-0.370587 & P3=-0.0115970 ;}
        nu_eff = 1.9986667e13 ;} // 15 micron
     end
     -3: begin
        P0=10.0287 & P1=17.8653 & P2=0.276804 &  P3=-0.0199558 ;}
        nu_eff = 2.4180000e17 ;} // 'effective' 1 keV
     end
     -4: begin
        P0=6.08087 & P1=10.8331 & P2=0.276802 &  P3=-0.0199597 ;}
        nu_eff = 1.2090000e18 ;} // 'effective' 5 keV
     end
     else: stop
  endcase

  if MODE LT 0 then lband = P0* 10.^(P3*x) + P1* 10.^(P2*x) ;
  return, 10.^(log_l_bol)/lband
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_qso_hrh07_lf_crosssec, nu
  sigma = 0.                    ;
  metallicity_over_solar = 1.   ;
  keV_in_Hz = 2.418e17          ;
  c_light = 2.998e8             ;
  micron  = 1.0e-6              ;

  ; For optical-IR regions, we use the Pei numerical approximations below.
  ;
  ; xsi = tau(lambda)/tau(B) is the ratio of extinction at lambda to the 
  ;    extinction in the B-band. 
  ; k = 10^21 (tau_B / NH)   (NH in cm^2) gives the dimensionless gas-to-dust
  ;    ratio, with k=0.78 for MW, k=0.16 for LMC, k=0.08 for SMC.
  ;    k is INDEPENDENT of the grain properties, and seems to scale rougly
  ;    linearly with metallicity
  ; so, for now, assume solar metallicity and k = k_MW = 0.78. we can rescale later.
  ;
  ; tau_B = ( NH / (10^21 cm^-2) ) * k --> SIGMA_B = k*10^-21  cm^2
  ; tau_lambda = xsi * tau_B --> SIGMA = xsi * SIGMA_B
  ;
  ; k = 0.78 for the MW
  ; k = 0.08 for the SMC, approximately in line with the MW/LMC/SMC metallicity 
  ;  sequence, so we take a k_MW then scaled by the metallicity
  k_dust_to_gas = 0.78 * metallicity_over_solar ;
  lambda_in_microns = c_light / nu / micron ;
  if  nu LT 0.03*keV_in_Hz then begin
     a = [185., 27., 0.005, 0.010, 0.012, 0.030] ;
     l = [0.042, 0.08, 0.22, 9.7, 18., 25.] ;
     b = [90., 5.50, -1.95, -1.95, -1.80, 0.00] ;
     n = [2.0, 4.0, 2.0, 2.0, 2.0, 2.0] ;
     R_V = 2.93          ;
     xsi = 0.
;     for i=0,5 do $
;        xsi += a[i] / ( (lambda_in_microns/l[i])^n[i] + (l[i]/lambda_in_microns)^n[i] + b[i] ) ;
     xsi = total( a / ( (lambda_in_microns/l)^n + (l/lambda_in_microns)^n + b) ) 
     sigma += xsi * k_dust_to_gas * 1.0e-21 ;
  endif

  ; For 0.03 keV < E < 10 keV  
  ;   (7.2e15 < nu[Hz] < 2.4e18  or   1.2 < lambda[Angstroms] < 413)
  ;   we use the photoelectric absorption cross sections of 
  ;   Morrison & McCammon (1983)
  ;     NOTE: these assume solar abundances and no ionization, 
  ;             the appropriate number probably scales linearly with both
  ;   (this is all for the COMPTON THIN regime)
  if  nu GT (0.03*keV_in_Hz*1.362/3.0) then begin ;// above Lyman edge
     stop
     sigma += morrison_photoeletric_absorption(nu/keV_in_Hz) ;
  endif

  ; Floor in cross-section set by non-relativistic (achromatic) Thompson scattering
  ;  (technically want to calculate self-consistently for the ionization state of the 
  ;   gas, but for reasonable values this agrees well with more detailed calculations 
  ;   including line effects from Matt, Pompilio, & La Franca; since we don't know the 
  ;   state of the gas (& have already calculated the inner reflection component), this
  ;   is the best guess)
  sigma += 6.65e-25 ;;

return, sigma;
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_qso_hrh07_lf, z, LOG_L_BOL=log_l_bol, MODE=mode, M_AB_grid=M_AB_grid

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'phi = x_qso_hrh07_lf(z) [v1.0]'
    return, -1
  endif 

  c = x_constants()

  ;; Lumniosties for evaluation (in L_sun)
  d_log_l_bol = 0.1
  log_L_bol = 7.0 + d_log_l_bol*dindgen(111L) + 0.008935 
  nlogL = n_elements(log_L_bol)
  

  ;; ;;;;;;;;;;;;;
  ;; Band luminosities
  if not keyword_set(MODE) then MODE = 0  ;; 0=Bol, -1=B-band
  case MODE of 
     0: print, 'x_qso_hrh07_lf: Returning Bolometric'
     -1: print, 'x_qso_hrh07_lf: Returning M_B'
     else: stop ;; not tested
  endcase

  l_band_grid = alog10(x_qso_hrh07_lf_lband(log_L_bol, MODE, NU_EFF=nu_eff))
  M_AB_grid = -2.5*l_band_grid + 2.5*alog10(nu_eff) - 32.38265724887536 ;
  ;// AB_nu for l_band_grid=nuLnu/L_sun_bol and nu in Hz


  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Full fit of HRH07
  P0=-4.8250643 
  P1=13.035753 
  P2=0.63150872  
  P3=-11.763560 
  P4=-14.249833                 
  P5=0.41698725 
  P6=-0.62298947 
  P7=2.1744386 
  P8=1.4599393 
  P9=-0.79280099                
  P10=0. 
  P11=0.
  P12=0.         
  P13=0.                
  P14=0.

  xsi = alog10((1. + z)/(1. + 2.)) 
  phi_star = P0          
  beta_min = 1.3         

  ;; // phi_star is constant as a function of redshift
  l_star   = P1 + P2*xsi + P3*xsi*xsi + P4*xsi*xsi*xsi ;
  ;; // l_star evolves (luminosity evolution) as a cubic polynomial in the 
  ;; //    convenient variable xsi=log_{10}((1+z)/(1+z_ref)) with z_ref=2
  gamma_1  = P5 * 10.^(xsi*P6) 
  ;; // gamma_1 the faint-end slope -- optionally evolves with redshift
  gamma_2  = 2.0 * P7 / (10.^(xsi*P8) + 10.^(xsi*P9)) ;
  ;; // gamma_2 the bright-end slope -- optionally evolves with redshift
  gamma_2 = (gamma_2 > beta_min) 
  ;; // cap the bright-end evolution to prevent unphysical divergence
  x = log_L_bol - l_star        ;
  log10_phi = phi_star - alog10(10.^(x*gamma_1) + 10.^(x*gamma_2)) 
  ;;  This gives the Bolometric Phi
  sav_bol_grid = 10.^log10_phi

  if MODE EQ 0 then return, sav_bol_grid  ;; Intrinsic bolometric

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; All that follows is to convert to "Observed" Phi

  ;; ;;;;;;;;;;;
  ;; Jacobian
  x = log_L_bol - 10.    ;
  lband = 0.             ;
  a=4 & b=5
  case mode of
     0:  begin
        P0=0. & P1 = 0 & P2 = 0 & P3 = 0
        jacob = 1.
     end
     -1: begin
        P0=8.99833 & P1=6.24800 & P2=-0.370587 & P3=-0.0115970
     end
     -2: begin
        P0=10.6615 & P1=7.40282 & P2=-0.370587 & P3=-0.0115970
     end
     -3: begin
        P0=10.0287 & P1=17.8653 & P2=0.276804 &  P3=-0.0199558
     end
     -4: begin
        P0=6.08087 & P1=10.8331 & P2=0.276802 &  P3=-0.0199597
     end
     else: stop
  endcase
  D1 = P0*(1.+P3)* 10.^(P3*x) + P1*(1.+P2)* 10.^(P2*x) ;
  D2 = P0* 10.^(P3*x) + P1* 10.^(P2*x) ;
  if mode LT 0 then jacob = D1/D2                 
  phi_bol_grid = 10.^log10_phi * jacob


  ;; ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Bolometric corrections
  ;double lb0,sig0,prefac,expfac ;
  phi_grid_out = dblarr(nlogL)
  ;sig0 = l_band_dispersion(l_bol_grid[i_lbol],nu) ;
  x = log_L_bol - 9.            ;
  case MODE of 
     0: sig0 = replicate(0.01, nlogL)
     -1: begin
        s0 = 0.08 & beta = -0.20 & s1 = 0.065 ;}
     end
     -2: begin 
        s0 = 0.03 & beta = -0.10 & s1 = 0.095 ;}
     end
     -3: begin 
        s0 = 0.01 & beta =  0.10 & s1 = 0.060 ;}
     end
     -4: begin 
        s0 = 0.04 & beta =  0.05 & s1 = 0.080 ;}
     end
     else: stop
  endcase
  if MODE LT 0 then sig0 = s0 * 10.^(beta*x) + s1

  for i_lbol=0,nlogL-1 do begin
     lb0  = l_band_grid[i_lbol] ;
     
     expfac = -0.5/(sig0[i_lbol]*sig0[i_lbol])  ;
     prefac = 1./(sig0[i_lbol] * sqrt(2.0*!pi)) * phi_bol_grid[i_lbol] * d_log_l_bol ;

     ;; Sum it
     phi_grid_out += prefac * exp(expfac*(lb0-l_band_grid)*(lb0-l_band_grid)) ;
  endfor

  phi_bol_grid = phi_grid_out
  phi_grid_out[*] = 0.
  ;printcol, l_band_grid, M_AB_grid, log_L_bol, phi_grid_out
  ;stop

  ;; ;;;;;;;;;;;;;;;
  ;; NH Column density correction
  ;; /* now, convolve over the distribution of column densities to get the 
  ;; //   post-obscuration bolometric QLF, adopting the observed distribution from 
  ;; //   Ueda et al. (2003)
  ;; // this considers a flat NH function in each of several bins
  ;; //		;; need to define a "lower limit" NH --> Ueda et al. use 20, our optical 
  ;; //		;;   calculations seem to suggest that's about right, so fine
  ;; */
  NH_MIN  = 20.0                ;
  NH_MAX  = 25.0                ;
  D_NH    = 0.01                ;
  N_NH    = long( (NH_MAX-NH_MIN)/D_NH + 1.0) ;
  NH = NH_MIN + dindgen(N_NH)*D_NH ;
  case mode of 
     0: tau = dblarr(N_NH)
     -1: tau = 10.^(NH) * x_qso_hrh07_lf_crosssec(2.998d8/(4400.0e-10)); // call at nu_B
     -2: tau = 10.^(NH) * x_qso_hrh07_lf_crosssec(2.998d8/(15.0e-6)) ;	 // call at 15microns
     else: stop
  endcase

  ;printcol, NH, tau
  ;stop

  ;; loop over the LF and attenuate everything appropriately
  n0 = 0                        ;
  eps = 1.7                     ;
  psi_44 = 0.47                 ;
  beta_L = 0.10                 ;
  psi_max = (1.+eps)/(3.+eps)   ;
  for iNH=0,N_NH-1 do begin
     NH_0 = NH[iNH] ;		
     ;; // need to interpolate to lay this over the grid already set up
     l_obs_grid = l_band_grid - tau[iNH]/alog(10.) ;
     phi_obs_grid = phi_bol_grid
     n0 = 0                     ;
     for i_lbol=0,nlogL-1 do begin
        li = l_band_grid[i_lbol] ;
        while((li GE l_obs_grid[n0+1]) AND ((n0+1) LT (nlogL-1))) do n0+=1 ;
        if (n0+1 LT nlogL)  then begin
           ;; // interpolate between the two l_obs_grid points
           p1 = alog10(phi_obs_grid[n0]) ;
           p2 = alog10(phi_obs_grid[n0+1]) ;
           p0 = p1 + (p2-p1) * ((li - l_obs_grid[n0])/(l_obs_grid[n0+1] - l_obs_grid[n0])) ;
        endif else begin
           ;; // extrapolate out
           p1 = alog10(phi_obs_grid[nlogL-2]) ;
           p2 = alog10(phi_obs_grid[nlogL-1]) ;
           p0 = p1 + (p2-p1) * ((li - $
                                 l_obs_grid[nlogL-2])/(l_obs_grid[nlogL-1] $
                                                             - l_obs_grid[nlogL-2])) ;
        endelse
        
        L_HX  = alog10(x_qso_hrh07_lf_lband(log_L_bol[i_lbol],-4)) ;
        psi = psi_44 - beta_L * (L_HX + alog10(3.9) + 33.0 - 44.0) ;
        if (psi LT 0.) then psi = 0.  ;
        if (psi GT psi_max) then psi = psi_max ;

        f_low = 2.0 - ((5.+2.*eps)/(1.+eps))*psi ;
        f_med = (1./(1.+eps))*psi ;
        f_hig = (eps/(1.+eps))*psi ;
        f_compton = f_hig       ;
        f_low = f_low / (1. + f_compton) ;
        f_med = f_med / (1. + f_compton) ;
        f_hig = f_hig / (1. + f_compton) ;	
        f_NH = 0.0              ;
        if ((NH_0 LE 20.5)) then f_NH = f_low ;
        if ((NH_0 GT 20.5) AND (NH_0 LE 23.0)) then f_NH = f_med ;
        if ((NH_0 GT 23.0) AND (NH_0 LE 24.0)) then f_NH = f_hig ;
        if ((NH_0 GT 24.0)) then f_NH = f_compton ;
        dN_NH = f_NH * D_NH     ;
        phi_grid_out[i_lbol] += 10.^p0 * dN_NH ;
     endfor
  endfor
  lsun=alog10(3.9)+33.           ;
  l_band_grid += lsun;
  log_L_bol  += lsun   ;
  ;printcol, l_band_grid, M_AB_grid, log_L_bol, phi_grid_out
  ;stop

  return, phi_grid_out
end
