;+
; NAME:
;   voigtwal
;
; PURPOSE:
;   Calculate a single voigt profile, given wavelength, absorber, and
;    an array of lines.
;
; CALLING SEQUENCE:
;   tau = voigtwal(wave,abs,lines)
;
; INPUTS:
;   wave       - Array of wavelengths to realize voigt profile
;   abs        - single absorber
;   lines      - line(s) structure to compute voigt profile, if mutliple,
;                   the entire list will be looped through.
;
; OPTIONAL INPUTS:
;   FWHM      - Resolution of instrument in PIXELS of input array
;
; OUTPUTS:
;   tau        - Optical depth at each wavelength
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;  Line is a structure with:
;
; ** Structure <81a5f94>, 4 tags, length=24, refs=3:
;    ION             STRING    'H I'
;    WAVE            DOUBLE           1215.6701
;    F               FLOAT          0.416400
;    GAMMA           FLOAT       6.26500e+08
;
;  Abs must have at least :
;
;** Structure ABS, 10 tags, length=44:
;   N               FLOAT           15.0000   	log 10 Column density
;   B               FLOAT           5.00000     velocity dispersion (km/s)
;   Z               FLOAT           2.00000	redshift
;
; EXAMPLES:
;
; BUGS: ;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   25-May-2000  Created by S. Burles, FNAL/IAP
;-
;------------------------------------------------------------------------------

function x_voigt, wave, lines, VELO=velo, FWHM=fwhm

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'tau (fx) = x_voigt( wave (vel), lines, /VELO, SIGMA= [v1.1]'
    return, -1
  endif 

  if keyword_set (VELO ) then stop  ;; Not ready for VELO right now
  if not keyword_set( FWHM ) then stop  ;; Need to deal with tau!
  

  npix = n_elements(wave)
  nsub = round((alog10(wave[npix-1]) - alog10(wave[0])) / 1.449E-6) + 1

  subwv = 10^(alog10(wave[0]) + dindgen(nsub)*1.449E-6)

  ;; Renormalize b
  bnorm = lines.b/299792.4581

  for i=0, n_elements(lines) - 1 do begin
      line = lines[i]
      vd = line.b/ (line.wrest * 1.0e-13)
      
      if not keyword_set( VELO ) then $
        vel = abs((subwv/(line.wrest*(1.0+line.zabs)) - 1.0)  / bnorm[i]) $
      else stop ; vel = wave / bnorm[i]  -- Not sure if this is valid anymore
      
      a = line.gamma / (12.56637 * vd)
      
      calc1 = where(vel GE 19.0) 
      calc2 = where(vel LT 19.0) 
      
      vo = vel*0.0	
      IF (calc1[0] NE -1) then begin 
          vel2 = vel[calc1]*vel[calc1]
          hh1 = 0.56419/vel2 + 0.846/(vel2*vel2)
          hh3 = -0.56 / (vel2 * vel2) 
          vo[calc1] = a * (hh1 + a * a * hh3) 
      ENDIF
      if (calc2[0] NE -1) then vo[calc2] = voigt(a,vel[calc2]) 
      
      thistau = 0.014971475*(10.0^line.N)*line.f*vo/vd
      if keyword_set(tau) then tau = tau + thistau $
      else tau = thistau
  endfor

  ;; Flux
  subfx = exp(-tau)

  ;; Create Gaussian Matrix
;  gmatr = dblarr(nsub, npix)
  ;; Create dwave
  dwv = wave - shift(wave,1)
  dwv[0] = dwv[1]

  ;; Final fx
  fx = fltarr(npix)

  soname = filepath('libxmath.so', $
                    root_dir=getenv('XIDL_DIR'), $
                    subdirectory='/lib')
  retval = call_external(soname, 'x_gsmooth', $
                         npix, nsub, double(wave), double(subwv), $
                         double(dwv), float(subfx), float(fx), $
                         float(FWHM), /UNLOAD)

  return, fx

end
	
        
