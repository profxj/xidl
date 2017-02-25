;+ 
; NAME:
; cldy_galx_sed
;   Version 1.1
;
; PURPOSE:
;   Given a starburst SED (wavelength, specific luminosity) scaled to
;   SFR, this code produces a Cloudy input spectrum
;
; CALLING SEQUENCE:
;   
;   cldy_galx_sed, sed_fil, outfil, SFR=, DISTANCE=, HM_z=
;
; INPUTS:
;   sed_fil  - FITS binary table with a structure tagged WAVE, STELLAR
;   distance - Distance from the galaxy in kpc  
;
; RETURNS:
;
; OUTPUTS:
;   outfil  - Cloudy output file
;
; OPTIONAL INPUTS:
;  SFR= --  Scaling of the luminosity;  SED is assumed scaled to
;          1Msun/yr
;  HM_z= -- Redshift to include a HM spectrum 
;  HM_SCALE=  -- Scale the HM spectrum by a factor
;  CUBA_FILE= -- CUBA output file 
;  FESC= -- Escape fraction at 912Ang (default=1)
;
; OPTIONAL OUTPUTS:
;  ENERGY=
;  LOGFNU=
;
; COMMENTS:
;
; EXAMPLES:
; cldy_galx_sed
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   31-Jul-2010 Written by JXP
;-
;------------------------------------------------------------------------------
pro cldy_galx_sed, sed_fil, distance, outfil, GEOM = GEOM $
                   , HM_z = HM_z, SFR = sfr $
                   , CUBA_FILE = cuba_file $
                   , HM_SCALE = hm_scale, ENERGY = energy, LOGFNU = logfnu $
                   , FESC = fesc $
                   , F912 = F912, HNUBAR = HNUBAR, CALCU = CALCU $
                   , PHI = PHI, E_MAX = E_MAX, GAMMA = GAMMA

  LN10 = alog(10.0D)
  c = x_constants()
  IF NOT KEYWORD_SET(E_MAX) THEN E_MAX = 1.0d10
  if  N_params() LT 3  then begin 
     print, 'Syntax - ' +$
            'cldy_galx_sed, sed_fil, distance (kpc), outfil,geom, HM_z=, ' + $
            'SFR=, HM_SCALE=, CUBA_FILE=, ENERGY=, LOGFNU=, FESC= [v1.1]'
     return
  endif 
  ;; If you are adding in HM spectrum must specify geometry. 
  IF KEYWORD_SET(HM_z) THEN BEGIN
     ;; For an optically thick illuminated slab, use 2.0*pi to compute
     ;; the phi factor. Otherwise use 4pi 
     CASE geom OF
        'thin': GEOM_FACT = 4.0d*!dpi
        'thick': GEOM_FACT = 2.0d*!dpi
        ELSE: message, 'geom must be either thick or thin'
     ENDCASE
  ENDIF
  if not keyword_set(HM_SCALE) then HM_SCALE = 1.  
  if not keyword_set(SFR) then sfr = 1.  ; Msun/yr
  if not keyword_set(FESC) then fesc = 1.
  if not keyword_set(NHI_FESC) then NHI_FESC = 1d20  ;; See Cantalupo 2010, MNRAS, 403, L16
  IF NOT KEYWORD_SET(CUBA_FILE) THEN $
     CUBA_FILE = getenv('XIDL_DIR') + '/Cloudy/cuba_q1g01_bkgthick.out'

  ;; Constants
  c = x_constants()

  ;; Read SED
  sed = xmrdfits(sed_fil, 1)
  
  ;; Convert to Energy units
  energy = c.RYDLAM/sed.wave ;; Ryd
  minE = min(energy, max=maxE)
  nen = n_elements(energy)

  flambda = 10.d^sed.stellar * SFR / (4*!pi*(distance*c.kpc)^2)
  fnu = flambda / c.c * sed.wave * (sed.wave * 1e-8)  ;; erg s^-1 cm^-2 Hz^-1

  ;; Allow for the escape fraction
  if FESC LT 1. then begin
     sigma_nu = x_photocross(1, 1, energy*c.Ryd/c.eV) 
     idx = where(energy GE 1.)
     fnu[idx] = fnu[idx] * (fesc + (1-fesc) * exp(-1*sigma_nu[idx]*NHI_fesc) )
  endif

  ;; Sort
  srt = sort(energy)
  energy = energy[srt]
  fnu = fnu[srt]
  log_energy = alog10(energy)
  mini = min(abs(energy-1.0D), k)
  F912 = fnu[k]

  ;; HM
  if keyword_set(HM_z) then begin
     cldy_read_cuba, cuba_file, HM_z, WAVE=HM_wave, JNU=HM_jnu
     ;; Scale and energy
     HM_jnu = HM_jnu * HM_SCALE
     HM_energy = c.RYDLAM / HM_wave
     ;; Sort
     HM_srt = sort(HM_energy)
     HM_energy=HM_energy[HM_srt]
     HM_jnu=HM_jnu[HM_srt]
     ;; Uniqe
     iuniq = uniq(HM_energy)
     HM_energy = HM_energy[iuniq] 
     HM_jnu = HM_jnu[iuniq] 

     ;; Report the relative value over 1-4 Ryd
     scl_idx = where(energy GE 1. and energy LT 4.)
     HM_val = interpol(HM_jnu, HM_energy, energy[scl_idx])
     sigma_nu = x_photocross(1, 1, energy[scl_idx]*c.Ryd/c.eV) 
     sigma_nu = sigma_nu / max(sigma_nu) 
     mean_scl = total(fnu[scl_idx] / (HM_val*GEOM_FACT) * sigma_nu) / $
               total(sigma_nu)
     print, 'Mean(Galaxy/EUVB) [1-4 Ryd] = ', mean_scl

     ;; Add into stellar 
     minHME = min(HM_energy, max=maxHME)
     add = where(energy GE minHME and energy LE maxHME, nadd) 
     if nadd GT 0 then begin
        HM_add = interpol(HM_jnu, HM_energy, energy[add])
        fnu[add] = fnu[add] + HM_add*GEOM_FACT  
     endif

     ;; Tack on other energy values as appropriate
     low_HME = where(HM_energy LT minE, nlow)
     if nlow NE 0 then begin
        energy = [HM_energy[low_HME], energy]
        fnu = [HM_jnu[low_HME]*GEOM_FACT, fnu]
     endif
     hi_HME = where(HM_energy GT maxE, nhi)
     if nhi NE 0 then begin
        energy = [energy,HM_energy[hi_HME]]
        fnu = [fnu,HM_jnu[hi_HME]*GEOM_FACT]
     endif
     
  endif
     
  logfnu = alog10(fnu)
  nu = energy*c.Ryd/c.h
  mn = min(abs(energy-1.1), imn)
  print, 'cldy_galx_sed: Normalize with:  f(nu)=', logfnu[imn], $
         ' [', energy[imn], ']'

  ;; Output
  close, /all
  openw, 2, outfil
  printf, 2, 'interpolate (0.00001 -30.0)'  ;; Ryd, log(Fnu)

  ;; Interpolate at each wavelength!
  FOR ii = 0L, nen-1L DO BEGIN
     printf, 2, 'continue ('+strtrim(energy[ii], 2)+' '+ $
             string(logfnu[ii], FORMAT = '(f10.6)')+')' 
  ENDFOR
  printf, 2, 'continue (7400000 -50.0)'
  close, 2

  IF arg_present(CALCU) OR arg_present(PHI) THEN BEGIN
     blue_energy = where(energy GE 1.0D AND energy LE E_MAX, nblue $
                         , COMPLEMENT = zero_energy, NCOMP = nzero)
     integrand = fnu/c.h
     integrand[zero_energy] = 0.0d
     phi = LN10*int_tabulated(log_energy, integrand, /double, /SORT)
     U =  phi/1.0d/c.c
     calcU = alog10(U)
  ENDIF

  IF arg_present(GAMMA) OR arg_present(HNUBAR) THEN BEGIN
     blue_energy = where(energy GE 1.0D, nblue $
                         , COMPLEMENT = zero_energy, NCOMP = nzero)
     sigma_nu = hcross_section(energy)
     integrand = fnu*sigma_nu/c.h*1.0d12
     integrand[zero_energy] = 0.0d
     gamma = LN10*int_tabulated(log_energy, integrand, /double, /SORT)
     ;; int 4pi*J_nu/(h*nu)*(h(nu-nu_LL))dnu
     integ_heat = fnu*sigma_nu*(energy-1.0d)*1.0d12*c.Ryd/c.h
     num = LN10*int_tabulated(log_energy, integ_heat, /double, /SORT)
     HNUBAR = num/gamma
  ENDIF
  
  RETURN
  
end
