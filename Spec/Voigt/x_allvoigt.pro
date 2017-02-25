;+
; NAME:
;   x_allvoigt
;  Version 1.1
;
; PURPOSE:
;   Calculate a single voigt profile, given wavelength, absorber, and
;    an array of lines.  Warning:  This program is not exact!
;
; CALLING SEQUENCE:
;   tau = x_allvoigt(wave,lines, SIGMA=, MNDV=)
;
; INPUTS:
;   wave       - Array of wavelengths to realize voigt profile
;   lines      - Line(s) structure to compute voigt profile, if mutliple,
;                   the entire list will be looped through.
;  RETURNS:
;   tau        - Optical depth at each wavelength
;
; OPTIONAL INPUTS:
;   SIGMA=     - Resolution of instrument in pixels
;   MNDV=      - Range in velocity to calculate Voigt over [defalut:
;                1000 km/s]
;
; OUTPUTS:
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
;                Modified by JXP
;-
;------------------------------------------------------------------------------

function x_allvoigt, wave, lines, SIGMA=sigma, MNDV=mndv

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'tau (fx) = x_allvoigt( wave, lines, SIGMA=, MNDV= ) [v1.1]'
    return, -1
  endif 
  
  ;; Optional keywords
  if not keyword_set(MNDV) then mndv = 1000.
  if not keyword_set(NPIX) then npix = n_elements(wave)

  tau = fltarr(npix)

  ;; Renormalize b
  spl = 2.997924581d5
  bnorm = lines.b/299792.4581

  for i=0, n_elements(lines) - 1 do begin
      ;; Grab line
      line = lines[i]
      vd = line.b/ (line.wrest * 1.0e-13)
      
      ;; Set dv
      if line.N GT 19.0 AND line.wrest GT 1000. then dv = 25000. $
      else dv = mndv

      ;; Find line center
      mn = min(abs(wave-line.wrest*(1.+line.zabs)), pcen)
      if pcen EQ 0 or pcen EQ (npix-1) then begin
         delv = abs(spl*(wave[1]-wave[0])/wave[0])
         pmn = (pcen - round(dv/delv)) > 0
         pmx = (pcen + round(dv/delv)) < (npix-1)
      endif else begin
         ;; Line-center is encompassed
         delv = spl*(wave[pcen+1]-wave[pcen])/wave[pcen]
         pmn = (pcen - round(dv/delv)) > 0
         pmx = (pcen + round(dv/delv)) < (npix-1)
      endelse
      vel = abs((wave[pmn:pmx]/(line.wrest*(1.0+line.zabs)) - 1.0)  / bnorm[i]) 
      
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
      
      thistau = 0.014971475*(10.0d^line.N)*line.f*vo/vd
      tau[pmn:pmx] = tau[pmn:pmx] + thistau 
  endfor

  ;; Return tau or fx
  if not keyword_set(SIGMA) then return, tau else begin
      kern = gauss_kernel(sigma)
      fx = exp(- convol(tau,kern) )
;      fx = exp(- (convol(tau,kern) < 80.) )
      return, fx
  endelse

end
	
        
