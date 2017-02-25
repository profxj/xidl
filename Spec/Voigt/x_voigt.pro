;+
; NAME:
;   x_voigt
;  Version 1.1
;
; PURPOSE:
;   Calculate a single voigt profile, given wavelength, absorber, and
;    an array of lines.  This can be a very expensive routine to run.
;
; CALLING SEQUENCE:
;   tau = x_voigt(wave, lines, FWHM=)
;
; INPUTS:
;   wave       - Array of wavelengths to realize voigt profile
;   lines      - line(s) structure to compute voigt profile. If mutliple,
;                   the entire list will be looped through.
; RETURNS:
;  Normalized flux array with Voigt profile super-posed
;
; OPTIONAL INPUTS:
;   FWHM      - Resolution of instrument in PIXELS of input array
;   COVERING  -- Covering fraction of the gas  [default = 1]
;
; OUTPUTS:
;   tau        - Optical depth at each wavelength
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;  x_gsmooth
;
;
; REVISION HISTORY:
;   25-May-2000  Created by S. Burles, FNAL/IAP
;   2003         Modified by JXP to give 'exact' answer
;   21-Jun-2007  Modified by KLC to opt out of 'exact' answer with /nosmooth
;-
;------------------------------------------------------------------------------

function x_voigt, wave, lines, VELO=velo, FWHM=fwhm, VERBOSE=verbose, $
                  nosmooth=nosmooth, TAU=tau, SUBWV=subwv, COVERING=covering

if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'tau (fx) = x_voigt( wave (vel), lines, /VELO, FWHM= [v1.1]'
    return, -1
endif 

if keyword_set (VELO ) then stop  ;; Not ready for VELO right now
if not keyword_set( FWHM ) and not keyword_set(NOSMOOTH) then stop  ;; Need to deal with tau!
if not keyword_set( COVERING ) then covering = 1.


npix = n_elements(wave)
if keyword_set(nosmooth) then begin
    nsub = npix
    subwv = wave
endif else begin
    nsub = round((alog10(wave[npix-1]) - alog10(wave[0])) / 1.449E-6) + 1
    subwv = 10^(alog10(wave[0]) + dindgen(nsub)*1.449E-6)
endelse

;; Renormalize b
bnorm = lines.b/299792.4581

for i=0, n_elements(lines) - 1 do begin
    if keyword_set(VERBOSE) AND (i mod 50 EQ 0) then $
      print, 'x_voigt: i = ', i
    line = lines[i]
    vd = line.b/ (line.wrest * 1.0e-13)
    
    if not keyword_set( VELO ) then $
      vel = abs((subwv/(line.wrest*(1.0+line.zabs)) - 1.0)  / bnorm[i]) $
    else stop ; vel = wave / bnorm[i]  -- Not sure if this is valid anymore
    
    a = line.gamma / (12.56637 * vd)
    
    calc1 = where(vel GE 19.0,complement=calc2) 
;      calc2 = where(vel LT 19.0) 
    
    vo = vel*0.0	
    IF (calc1[0] NE -1) then begin 
        vel2 = vel[calc1]*vel[calc1]
        hh1 = 0.56419/vel2 + 0.846/(vel2*vel2)
        hh3 = -0.56 / (vel2 * vel2) 
        vo[calc1] = a * (hh1 + a * a * hh3) 
    ENDIF
    if (calc2[0] NE -1) then begin
      ;stop
      vo[calc2] = voigt(a,vel[calc2]) 
    endif
    
    thistau = 0.014971475*(10.0^line.N)*line.f*vo/vd
    if i GT 0 then tau = tau + thistau  else tau = thistau
endfor

;; Flux
subfx = exp(-tau)

;; Covering fraction
if covering LT 1. then begin
   subfx = 1. - COVERING*(1-subfx)
endif

;; Create Gaussian Matrix
;  gmatr = dblarr(nsub, npix)
;; Create dwave
dwv = wave - shift(wave,1)
dwv[0] = dwv[1]

;; Final fx
if keyword_set(nosmooth) then fx = subfx $
else begin
    fx = fltarr(npix)

    if keyword_set(VERBOSE) then print, 'x_voigt: Entering x_gsmooth'
    soname = filepath('libxmath.' + idlutils_so_ext(), $
                      root_dir=getenv('XIDL_DIR'),  $
                      subdirectory='/lib')
    retval = call_external(soname, 'x_gsmooth', $
                           npix, nsub, double(wave), double(subwv), $
                           double(dwv), float(subfx), float(fx), $
                           float(FWHM))
;                         float(FWHM), /UNLOAD)
endelse 

return, fx

end
	
        
