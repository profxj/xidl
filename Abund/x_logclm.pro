;+ 
; NAME:
; x_logclm
;   Version 1.1
;
; PURPOSE:
;    Convert linear column density to log values
;
; CALLING SEQUENCE:
;   x_logclm, ntot, sig, logN, logS
;
; INPUTS:
;   ntot - Linear column density
;   sig  - Error
;
; RETURNS:
;
; OUTPUTS:
;  logN - Log10 N
;  logS - Error in Log10(N)
;
; OPTIONAL KEYWORDS:
;  /REVERSE  -- Convert log values to linear
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
;   05-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_logclm, ntot, sig, logN, logS, REVERSE=reverse

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_logclm, ntot, sig, logN, logS, /REVERSE [v1.1]'
    return
  endif 

  ;; Optional Keywords
  if not keyword_set(REVERSE) then begin
      logN = alog10(ntot)
      lgvar = ((1.0 / (alog(10.0)*ntot))^2)*sig^2
      logS = sqrt(lgvar)
  endif else begin
      ntot = 10.^logN
      sig = logS * alog(10.) * Ntot
  endelse

  return
end
