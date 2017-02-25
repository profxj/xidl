;+ 
; NAME:
; x_calctb   
;   Version 1.1
;
; PURPOSE:
;    Calculate the temperature of a gas given its Doppler parameter
;
; CALLING SEQUENCE:
;   x_calctb, b, T, /LOGT, MA=
;
; INPUTS:
;  b -- Doppler parameter (km/s)
;
; RETURNS:
;
; OUTPUTS:
;  T -- Temperature (K)
;
; OPTIONAL KEYWORDS:
;  /LOGT -- Return logarithmic temperature
;  MA -- Mass of nucleus [1 for Hydrogen]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  x_calctb, 25., T, /LOGT
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_calctb, b, T, MA=MA, LOGT=logt

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_calctb, b, t, MA=, /LOGT, [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( MA ) then MA = 1.

;; Constants
;  c = x_constants(/cgs)
  c = x_constants()

  T = (b*1e5)^2 * MA * c.mp / c.k / 2.
  if keyword_set( LOGT ) then T = alog10(T)

  return
;    
end
