;+ 
; NAME:
; x_calctau0
;   Version 1.1
;
; PURPOSE:
;    Calculate the peak optical depth (tau_0) of a transition given its rest
;    wavelength, Doppler parameter and column desity.
;
; CALLING SEQUENCE:
;    tau0 = x_calctau0(wrest, b, N /REVERSE) [v1.0]'
;
; INPUTS:
;   wrest -- Wavelength (Ang)
;   b  - Doppler (km/s)
;   N  - Log column (or tau without log if /reverse)
;
; RETURNS:
;
; OUTPUTS:
;  tau0 - Peak optical depth (or log N for /Reverse)
;
; OPTIONAL KEYWORDS:
;  /REVERSE -- Return N instead of tau0 (assumes input of tau_0)
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
;   21-Apr-2006 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_calctau0, wrest, b, N, REVERSE=reverse

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'tau0 = x_calctau0(wrest, b, N /REVERSE) [v1.1]'
    return, -1
  endif 

  c = x_constants()

  ;; Constant
  getfnam, wrest, fval
  cnst = sqrt(!pi) * c.e^2 * (wrest*1d-8) * fval / (c.me*c.c*b*1e5)
  
  ;; 
  if keyword_set(REVERSE) then begin
      ;; Return N
      tau0 = alog10(float(N)/cnst)
  endif else begin
      ;; Return tau0
      tau0 = cnst * 10^double(N)
  endelse

  return, tau0
end
