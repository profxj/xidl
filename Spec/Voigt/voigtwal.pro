;+
; NAME:
;   voigtwal
;
; PURPOSE:
;   Calculate a single voigt profile, given wavelength, absorber, and
;    an array of lines.  I recommend x_voigt
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
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   25-May-2000  Created by S. Burles, FNAL/IAP
;-
;------------------------------------------------------------------------------

function voigtwal,wave,abs,lines

	bnorm = abs.b/299792.4581

        for i=0, n_elements(lines) - 1 do begin
          line = lines[i]
 	  vd = abs.b/ (line.wave * 1.0e-13)

	  vel = abs((wave/(line.wave*(1.0+abs.z)) - 1.0)  / bnorm)
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
 
          thistau = 0.014971475*(10.0^abs.n)*line.f*vo/vd
	  if keyword_set(tau) then tau = tau + thistau $
	  else tau = thistau
        endfor
      return, tau
end
	
        
