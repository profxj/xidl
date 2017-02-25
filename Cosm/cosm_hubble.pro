;+ 
; NAME:
; cosm_hubble
;
; PURPOSE:
;   Calculates Hubbles constant at arbitrary redshift
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;   z  -- Redshift
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /INIT  -- Initializes the cosmology to the default values
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hubb = cosm_hubble(2., /INIT)
;   
;
; PROCEDURES CALLED:
;   cosm_common
;
; REVISION HISTORY:
;   22-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

function cosm_hubble, z, H0=h0, INIT=init, _EXTRA=extra

  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, cosm_L, cosm_r, sigma_8

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'hubb = cosm_hubble(z, H0=, /INIT) [v1.1]'
    return, -1
  endif 

  if keyword_set( INIT ) or keyword_set(EXTRA) then $
    cosm_common, H0=h0, _EXTRA=extra
  if keyword_set( H0 ) then cosm_h = h0

  cosm_hubble = sqrt(cosm_L + (cosm_K)*(1+z)^2 + $
                     cosm_dm*(1+z)^3 + cosm_r*(1+z)^4)

  return, cosm_hubble*cosm_h

end

