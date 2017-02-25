;+ 
; NAME:
; x_qsonum
;   Version 1.1
;
; PURPOSE:
;    Returns number of qso's per sq deg for a redshift and magnitude
;    interval.   For now assumes 2DF and b_J magnitudes
;
; CALLING SEQUENCE:
;   
; INPUTS:
;  mag_lim -- Limiting magnitude of the survey
;  z_lim   -- Redshift interval for the QSOs
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
;   07-Jan-2006 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_qsonum, mag_lim, z_lim 

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'nqso = x_qsonum(mag_lim, z_lim) [v1.0]'
    return,-1
  endif 
  
  if z_lim[0] LT 1.8 then stop
  if z_lim[1] GT 3.5 then begin
      print, 'Warning: Am assuming no QSO evolution at z>3.5!'
      stop
  endif

  ;; 2DF number counts
  ;; z=1.8 to 2.1
  mag = -24.25 - 0.5*findgen(9)  ;; bins of 0.5 mag
  nqso = [26, 897., 898., 618, 333, 120, 43, 9, 2]  ; I bumped 2 -> 9
  effarea = 673.4  ; sq deg

  ;; 

  return
end
