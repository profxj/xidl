;+ 
; NAME:
; cldy_calcu
;   Version 1.0
;
; PURPOSE:
;    Creates a Uplot for a given NHI, FeH, nH
;
; CALLING SEQUENCE:
;   
; cldy_Uplot, grid, NHI, FeH, nH, ions
;
; INPUTS:
;   grid  - CLOUDY grid
;   NHI - Can be an array of values
;   FeH
;   nH
;   ions  - Array of [Z,ion] vectors
;
; RETURNS:
;   
;
; OUTPUTS:
;   Creates a Plot
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   cldy_calcu, grid, 19.0d, -1.0d, -1.0d, [[14,2], [14,3]]
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cldy_calcu, grid, z, nH, logU, NRM=nrm

;
  if  N_params() LT 3  then begin 
      print, 'Syntax - ' +$
        'cldy_calcu, grid, z, nH, [logU], NRM= [v1.0]'
      return
  endif 

; Optional keywords

  if not keyword_set( NRM ) then nrm = 1.e-23

  ;; Get redshifts
  close, /all
  openr, 1, grid

  dumc = ' '
  for i=0L,31 do readf, 1, dumc
  zval = fltarr(10)
  readf, 1, zval, FORMAT='(10x,10f10.4)' 
  close, /all
  

  ;; Grep flux
  readcol, grid, wave, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, skipline=34
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

  Phi = total(4*!pi*NRM * nflux[0:imn]*dnu[0:imn] / (6.63e-27*nu[0:imn]))

  ;; U
  logU = alog10(Phi / (nH * 3e10))

end

