;+ 
; NAME: 
; ew_teff_lyman   
;    Version 1.1
;
; PURPOSE:
;    Calculate tau_effective for the Lyman series using the EW
;    approximation (e.g. Zuo 93)
;
; CALLING SEQUENCE:
;   
; teff = ew_teff_lyman(wave,zem) 
;
; INPUTS:
;  ilambda -- Observed wavelength (scalar)
;  zem -- Emission redshift of the source [sets which Lyman lines are
;         included]
;
; RETURNS:
;  teff -- Total effective opacity of all lines contributing
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  bva=    -- Characteristics Doppler parameter for the Lya forest
;             [Options: 24, 35 km/s]
;  NHI_MIN -- Minimum log HI column for integration [default = 11.5]
;  NHI_MAX -- Maximum log HI column for integration [default = 22.0]
;  /FNZ -- f(N) is in dz not dX
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   teff = ew_teff_lyman(3400., 2.4)
;
; PROCEDURES/FUNCTIONS CALLED:
;  
;
; REVISION HISTORY:
;   May-2011 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;
;;;;
function ew_teff_lyman, ilambda, zem, EW_SPLINE=EW_SPLINE, bval=bval, $
                        LYMAN=lyman, POWERFN=powerfn, N_eval=n_eval, $
                        TEFF_LYMAN=teff_lyman, FNZ=fnz, DEBUG=debug, $
                        NHI_MIN=nhi_min, NHI_MAX=nhi_max, WREST=wrest, $
                        NLYMAN=nlyman, INTGRND=intgrnd, LGNVAL=lgnval

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'teff = ew_teff_lyman(wave, zem, /DEBUG, NHI_MIN=, NHI_MAX= [v1.0])'
    	return, -1
  endif 

  if not keyword_set(bval) then bval = 35. ; km/s
  if not keyword_set(N_eval) then N_eval = 5000L
  if not keyword_set(NHI_MIN) then NHI_MIN = 11.5
  if not keyword_set(NHI_MAX) then NHI_MAX = 22.0 

  ;; Parse input array
  if n_elements(ilambda) GT 1 then stop; return, 0
  lambda = ilambda[0]

  ;; Read in EW spline (if needed)
  if not keyword_set(EW_SPLINE) then begin
     case round(bval) of 
        24: EW_FIL = getenv('XIDL_DIR')+'/IGM/EW_SPLINE_b24.fits'
        35: EW_FIL = getenv('XIDL_DIR')+'/IGM/EW_SPLINE_b35.fits'
        else: stop
     endcase
     EW_SPLINE = xmrdfits(EW_FIL, 1, /silent)
  endif
     
  ;; Assumed Line list
  wrest = [1215.6701d,  1025.7223,       972.53680,       949.74310,       937.80350, $
           930.74830,  926.22570,       923.15040,       920.96310,       919.35140, $
           918.12940,  917.18060,       916.42900,       915.82400,       915.32900, $
           914.91900,  914.57600,       914.28600,       914.03900,       913.82600, $
           913.64100,  913.48000,       913.33900,       913.21500,       913.10400, $
           913.00600,  912.91800,       912.83900,       912.76800,       912.70300, $
           912.64500]

  ;; Initialize f(N)
  if not keyword_set(POWERFN) then init_powerfn, powerfn

  ;; Find the lines
  gd_Lyman = where(lambda/(1+zem) LT wrest, nlyman)
  if nlyman EQ 0 then return, -1

  
  ;; N_HI grid
  lgNval = NHI_MIN + (NHI_MAX-NHI_MIN)*dindgen(N_eval)/(N_eval-1) ;; Base 10 
  dlgN = lgNval[1]-lgNval[0]
  Nval = 10.d^lgNval
  teff_lyman = dblarr(nlyman)

  for qq=0L,nlyman-1 do begin  ;; Would be great to do this in parallel... 
                               ;;; (Can pack together and should)
     ;; Redshift
     zeval = (lambda / wrest[qq]) - 1
     if not keyword_set(FNZ) then $
        dxdz = abs(cosm_xz(zeval-0.1, /silent)- $
                   cosm_xz(zeval+0.1,/silent)) / 0.2 $
        ;dxdz = abs(cosm_xz(zeval-0.1, /W05MAP,/silent)- $
        ;           cosm_xz(zeval+0.1,/W05MAP,/silent)) / 0.2 $
     else dxdz = 1. ;; Code is using f(N,z)

     ;; Get EW values (could pack these all together)
     if abs(EW_SPLINE[qq].wrest-wrest[qq]) GT 1e-3 then stop
     restEW = spl_interp(EW_SPLINE[qq].NHI, EW_SPLINE[qq].EW, EW_SPLINE[qq].splint, lgNval) ;; Ang

     ;; dz
     dz = restEW * (1+zeval) / wrest[qq]

     ;; Evaluate f(N,X)
     log_fnX = eval_powerfn(POWERFN, lgNval, zeval)

     ;; Sum
     intgrnd = 10.d^(log_fnX) * dxdz * dz * Nval
     teff_lyman[qq] = total(intgrnd) * dlgN * alog(10.)
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
     
  endfor
  
  junk = check_math()

  return, total(teff_lyman)
end
