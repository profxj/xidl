;+ 
; NAME:
; voigtq
;    Version 1.0
;
; PURPOSE:
;    Creates a quick Voigt profile using code written by SB.
;   Not recommended for usage
;
; CALLING SEQUENCE:
;   
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------
;  

;  Line is a structure with:
;
; ** Structure <81a5f94>, 4 tags, length=24, refs=3:
;    ION             STRING    'H I'
;    WAVE            DOUBLE           1215.6701
;    F               FLOAT          0.416400
;    GAMMA           FLOAT       6.26500e+08
;
;  Abs is a structure with:
;
;** Structure ABS, 10 tags, length=44:
;   ION             STRING    'H I'
;   NUM             LONG                -1
;   N               FLOAT           15.0000   	log 10 Column density
;   N_ERR           FLOAT           0.00000
;   B               FLOAT           5.00000     velocity dispersion (km/s)
;   B_ERR           FLOAT           0.00000
;   Z               FLOAT           2.00000	redshift
;   Z_ERR           FLOAT           0.00000
;   NB              FLOAT          -2.00000     N-B correlation coefficient
;   REGIONS         POINTER   <NullPointer>
;

function voigtq,wave, abs, line

	bnorm = abs.b/299792.4581
	vd = abs.b/ (line.wave * 1.0e-13)

	vel = abs((wave/(line.wave*(1.0+abs.z)) - 1.0)  / bnorm)
        a = line.gamma / (12.56637 * vd)

	calc1 = where(vel GE 10.0) 
	calc2 = where(vel LT 10.0) 

	vo = vel*0.0	
	IF (calc1[0] NE -1) then begin 
	  vel2 = vel[calc1]*vel[calc1]
	  hh1 = 0.56419/vel2 + 0.846/(vel2*vel2)
	  hh3 = -0.56 / (vel2 * vel2) 
	  vo[calc1] = a * (hh1 + a * a * hh3) 
        ENDIF
	if (calc2[0] NE -1) then vo[calc2] = voigt(a,vel[calc2]) 
  
	tau = 0.014971475*(10.0^abs.n)*line.f*vo/vd
	return, exp(-tau)
end
	
        
