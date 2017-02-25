;+ 
; NAME:
; x_relvel
;  Version 1.1
;
; PURPOSE:
;  Calculate redshifts [default] or velocties using special rel
;  All velocities are km/s and positive means lower redshift
;
; CALLING SEQUENCE:
;   [z2 or v] = x_relvel( z1, [v or z2], /REVERSE)
;
; INPUTS:
;    z1 -- Redshift of 'rest'
;
; RETURNS:
;    v or z2 -- z2 is default given a velocity (km/s)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   
;
; OPTIONAL OUTPUTS:
;  AMIN -- Mpc per arcmin at z
;
; COMMENTS:
;
; EXAMPLES:
;   print, x_relvel(3., 3000.)  [Answer = 2.9601706]
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Aug-2006 Written by JXP
;-
;------------------------------------------------------------------------------

function x_relvel, z1, v, REVERSE=reverse

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             '[z2 or v] = x_relvel(z1, v or z2), /REVERSE [v1.0]'
    return, -1
  endif 

  ;; Standard
  c = x_constants()
  if not keyword_set(REVERSE) then begin
      ;; b
      b = (v*1e5)/c.c
      ;; R
      R = sqrt((1+b)/(1-b))
      ;; Finally
      z2 = (1+z1)/R - 1
  endif else begin
      ;; Velocity
      R = (1.+z1) / (1.+v)
      z2 = c.c * (R^2 - 1)/(1+R^2) /1e5 ; km/s
  endelse
      
  return, z2

end

