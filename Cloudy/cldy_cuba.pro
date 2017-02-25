;+ 
; NAME:
; cldy_cuba
;   Version 1.1
;
; PURPOSE:
;   Creates a Cloudy input file from a CUBA output file given a
;    redshift 
;
; CALLING SEQUENCE:
;   
;   cldy_cuba, fil, z, outfil, FIXG=fixg
;
; INPUTS:
;   fil  - CUBA output file
;   z    - Redshift
;
; RETURNS:
;
; OUTPUTS:
;   outfil  - Cloudy output file
;
; OPTIONAL INPUTS:
;  /FIXG -- I do not remember what this is for!
;  CALCU= -- Calculate the U parameter assuming n_H = 1 cm^-3
;
; OPTIONAL OUTPUTS:
; STRCT -- Stucture containing the wavelength and flux values
;
; COMMENTS:
;
; EXAMPLES:
; cldy_cuba, 0.35, '/u/xavier/Cloudy/Spec/Output/q1g0_z035.spec'
; CUBA_FILE='/u/xavier/Cloudy/Spec/Data/CUBA/Q1G0/bkgthick.out',
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  cldy_read_cuba
;
; REVISION HISTORY:
;   06-Nov-2003 Written by JXP
;   31-Jul-2010 Grabbed JFHs version
;   15-Nov-2016 Until can get JFHs version, make this work; KLC
;-
;------------------------------------------------------------------------------
PRO cldy_cuba, z, outfil $
               , CALCU = calcu, PHI = phi $
               , GAMMA = GAMMA $
               , J912 = J912, JNU_SCALE = JNU_SCALE $
               , CUBA_FILE=cuba_file, skip_stop=skip_stop, _extra=extra

  if not keyword_set(skip_stop) then $
     ;; For those of us lacking cldy_cuba_jfh.pro
     stop,'cldy_cuba stop: WARNING! Use cldy_cuba_jfh until further notice [JXP: 12-12-2010]'

  RYDLAM = 911.26705 ;; converts Rydbergs to Angstroms
  IF NOT KEYWORD_SET(CUBA_FILE) THEN $
     CUBA_FILE = getenv('XIDL_DIR') + '/Cloudy/cuba_q1g01_bkgthick.out'

  IF NOT KEYWORD_SET(JNU_SCALE) THEN JNU_SCALE = 1.0D
  
  cldy_read_cuba, cuba_file, z, WAVE = WAVE, JNU = JNU1, _extra=extra ; _extra includes /hm2012
  ;; re-scale the JNU
  JNU = JNU_SCALE[0]*JNU1
  ;; Interpolate at each wavelength!
  energy = RYDLAM/wave
  isort = sort(energy)
  energy = energy[isort]
  logjnu = alog10(JNU[isort])
  iuniq = uniq(energy)
  energy = energy[iuniq] 
  logjnu = logjnu[iuniq] > (-30.0)
  nen = n_elements(energy)
  
  openw, 2, outfil
  printf, 2, 'interpolate (0.00000001 -30.0)' ; span CLOUDY range exactly

  FOR ii = 0L, nen-1L DO BEGIN
     printf, 2, 'continue ('+strtrim(energy[ii], 2)+' '+ $
             string(logjnu[ii], FORMAT = '(f10.6)')+')' 
  ENDFOR
  printf, 2, 'continue (7400000 -30.0)'
  close, 2
 
  IF arg_present(CALCU) OR arg_present(PHI) THEN BEGIN
     ;; Sum up phi
     c = x_constants()
     nu = (2.9979246d18/wave)   ; (A/s)
     blue_wave = where(wave LE 912.0D, nblue)
     nu_blue = nu[blue_wave]
     jnu_blue = jnu[blue_wave]
     ;; restrict to unique values
     uni_ind = uniq(nu_blue, sort(nu_blue))
     nu_blue = nu_blue[uni_ind]
     jnu_blue = jnu_blue[uni_ind]
     lognu_blue = alog10(nu_blue)
     mini = min(abs(wave-912.0D), k)
     J912 = jnu[k]
     LN10 = alog(10.0D)
     integrand = 4*!dpi*jnu_blue/c.h
     phi = LN10*int_tabulated(lognu_blue, integrand, /double, /SORT)
     U =  phi/1.0d/c.c
     calcU = alog10(U)
  ENDIF

  IF arg_present(GAMMA) THEN BEGIN
     ;; Sum up phi
     c = x_constants()
     nu = (2.9979246d18/wave)   ; (A/s)
     blue_wave = where(wave LE 912.0D, nblue)
     nu_blue = nu[blue_wave]
     jnu_blue = jnu[blue_wave]
     ;; restrict to unique values
     uni_ind = uniq(nu_blue, sort(nu_blue))
     nu_blue = nu_blue[uni_ind]
     jnu_blue = jnu_blue[uni_ind]
     lognu_blue = alog10(nu_blue)
     mini = min(abs(wave-912.0D), k)
     J912 = jnu[k]
     LN10 = alog(10.0D)
     sigma_nu = hcross_section(lognu_blue)
     ;; output is in units of 1d12
     integrand = 4*!dpi*jnu_blue*sigma_nu/c.h*1.0d12
     gamma = LN10*int_tabulated(lognu_blue, integrand, /double, /SORT)
  ENDIF
 
  RETURN
END
