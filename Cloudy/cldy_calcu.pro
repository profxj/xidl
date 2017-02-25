;+ 
; NAME:
; cldy_calcu
;   Version 1.0
;
; PURPOSE:
;    Calculate the ionization parameter 'U' for a given redshift
;    an nH value assuming the Haardt & Madau (1996) spectrum
;
; CALLING SEQUENCE:
;  cldy_calcu, hm_fil, z, nH, logU, NRM=nrm
;
; INPUTS:
;   hm_fil - Haardt & Madau file of fluxes 
;   z - Redshif
;   nH - Hydrogen volume density [linear]
;
; RETURNS:
;   
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;  logU -- The log of the ionization parameter
;
; OPTIONAL KEYWORDS:
;   [NRM] - Normalization factor for H&M data  [default: 1e-23]
;
; COMMENTS:
;
; EXAMPLES:
;   cldy_calcu, hm_fil, z, nH, logU, NRM=
;
; PROCEDURES/FUNCTIONS CALLED:
;  readcol
;
; REVISION HISTORY:
;   02-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cldy_calcu, hm_fil, z, nH, logU, NRM=nrm

;
  if  N_params() LT 3  then begin 
      print, 'Syntax - ' +$
        'cldy_calcu, hm_fil, z, nH, [logU], NRM= [v1.0]'
      return
  endif 

; Optional keywords

  if not keyword_set( NRM ) then nrm = 1.e-23

  ;; Get redshifts
  close, /all
  openr, 1, hm_fil

  dumc = ' '
  for i=0L,31 do readf, 1, dumc
  zval = fltarr(10)
  readf, 1, zval, FORMAT='(10x,10f10.4)' 
  close, /all
  

  ;; Grep flux
  readcol, hm_fil, wave, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, skipline=34
  tflux = dblarr(10,167)
  tflux[0,*] = f1
  tflux[1,*] = f2
  tflux[2,*] = f3
  tflux[3,*] = f4
  tflux[4,*] = f5
  tflux[5,*] = f6
  tflux[6,*] = f7
  tflux[7,*] = f8
  tflux[8,*] = f9
  tflux[9,*] = f10
  
  ;; Interpolate to get flux(z)

  a = where( z GE zval and z lt shift(zval,1), na)
  if na NE 1 then stop
  iz = a[0]

  wgt = (z - zval[iz]) / (zval[iz+1]-zval[iz])
  
  flux = tflux[iz,*]*(1.-wgt) + wgt*tflux[iz+1,*]


  ;; Calculate Phi
  nu = 2.99e18 / wave ; (A/s)

  mn = min( abs(wave-912.), imn)
  nflux = flux / flux[imn]
  dnu = nu - shift(nu,-1)

  Phi = total(4*!pi*NRM * nflux[0:imn]*dnu[0:imn] / (6.63d-27*nu[0:imn]))

  ;; U
  logU = alog10(Phi / (nH * 3e10))

end

