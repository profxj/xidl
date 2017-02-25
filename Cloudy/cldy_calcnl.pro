;+ 
; NAME:
; cldy_calcnl
;   Version 1.1
;
; PURPOSE:
;   Creates a Cloudy input file from a CUBA output file given a
;    redshift 
;
; CALLING SEQUENCE:
;   
;   cldy_calcnl, fil, z, logU, Jnu, NHVAL=, NHH=, LVAL=
;
; INPUTS:
;   fil  - CUBA output file
;   z    - Redshift
;   logU    - Ionization parameter
;   J912  - Intensity at Lyman limit
;
; RETURNS:
;
; OUTPUTS:
;   NHV=   -  Value of the volume density
;   LVAL= -  Length of the absorber  (kpc)
;
; OPTIONAL INPUTS:
;   NHH=  -  NH value of the sightline (required for lval)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; cldy_calcnl, '/u/xavier/Cloudy/Spec/Data/CUBA/Q1G0/bkgthick.out',
; 0.5, -1.1, 6e-23, NHVAL=nhval
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   30-Wed-2004 Written by JXP
;-
;------------------------------------------------------------------------------
pro cldy_calcnl, fil, z, logU, J912, NHV=nhv, NHH=nhH, LVAL=lval

  if  N_params() LT 4  then begin 
      print, 'Syntax - ' +$
        'cldy_calcnl, fil, z, logU, J912, NHV=, NHH=, LVAL= [v1.1]'
      return
  endif 

  ;; Open Madau file
  close, /all
  openr, 1, fil
  zval = fltarr(10)
  wv = dblarr(432L)
  flux = dblarr(432L,10)
  dumf = dblarr(11)

  ;; Loop
  for ii=0L,999 do begin
      readf, 1, zval, FORMAT='(11x,10f11.4)' 
      for jj=0L,431 do begin
          readf, 1, dumf;, FORMAT='(10f11.4)'
          flux[jj,*] = dumf[1:*]
          wv[jj] = dumf[0]
      endfor
      if (z GE zval[0] OR z LE zval[9]) then break
  endfor

  close, 1

  ;; Spline at each wavelength!
  esv = 0.
  energy = dblarr(431L)
  jnu = dblarr(431L)
  cnt = 0L
  svwv = dblarr(431L)
  for ii=431L,0,-1 do begin

      egy = 912./wv[ii]
      if egy NE esv then begin
          jnu[cnt] = interpol(flux[ii,*], zval, z, /spline)
          energy[cnt] = egy
          svwv[cnt] = wv[ii]
          cnt = cnt + 1
      endif
;      if fx GT 1E-30 then jnu = alog10(fx) else jnu = -30.
;      if keyword_set(FIXG) and wv[ii] LT 50 then jnu = -30.
;      if energy NE esv then $  ; Multiple energies in a few spots
;        printf, 2, 'continue ('+strtrim(energy,2)+' '+ $
;        string(jnu,FORMAT='(f7.3)')+')'
      esv = egy
  endfor
  svwv = svwv[0:cnt-1]
  energy = energy[0:cnt-1]
  jnu = jnu[0:cnt-1]

  ;; Sum up phi
  nu = 2.99e18 / svwv ; (A/s)
  gdwv = where(svwv LT 912., ngd)
  mn = min( abs(svwv-912.), imn)
  njnu = jnu[gdwv] / jnu[imn]
  dnu = shift(nu[gdwv],-1) - nu[gdwv]
  dnu[ngd-1] = dnu[ngd-2]
  Phi = total(4*!pi*J912 * njnu*dnu / (6.63e-27*nu[gdwv]))

  ;; calculate nH
  nHV = Phi / (10^logU * 3e10)
  if keyword_set( NHH ) and arg_present( LVAL ) then begin
      cnst = x_constants()
      lval = 10^NHH / nHV / cnst.kpc  ; kpc
  endif
  
  return
  
end
