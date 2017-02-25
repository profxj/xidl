;+ 
; NAME:
; cosm_xz
;
; PURPOSE:
;    Calculate the cosmological distance X from z=0 to z=z
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;   z -- Redshift
;
; RETURNS:
;   x -- Cosmological pathlength
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   H0 -- Hubbles constant (km/s/Mpc)
;   OM -- Omega Dark Matter
;   OV -- Lambda
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
; PROCEDURES CALLED:
;  cosm_common
;  cosm_intxz
;  qromb
;
; REVISION HISTORY:
;   11-March-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function cosm_intxz, z
  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, $
     cosm_L, cosm_r, sigma_8

  intxz = (1.+z)^2 / sqrt(cosm_L + (cosm_K)*(1+z)^2 + cosm_dm*(1+z)^3)
  return, intxz
end
  

function cosm_xz, z, zmin=zmin, H0=h0, OM=om, OV=ov, NOINIT=noinit, $
                  exact=exact, _EXTRA=extra

  common cosmolgy_cmmn

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'xz = cosm_xz(z, zmin= H0=, OM=, OV=, /W05MAP, /NOINIT) [v1.1]'
    return, -1
  endif 
  
  if not keyword_set(zmin) then zmin = 0.

  if not keyword_set(NOINIT) and not keyword_set(cosm_h) then $
     cosm_common, H0=h0, Omegavac=OV, OmegaDM=OM, _EXTRA=extra

  if keyword_set(exact) and abs(cosm_dm+cosm_L-1.) lt 1.e-5 then begin 
     xz = 2.d/(3*cosm_dm) * (sqrt(cosm_dm*(1+z)^3 + cosm_L) - $
                             sqrt(cosm_dm*(1+zmin)^3 + cosm_L)) ;; Checked with Wolfram
  endif else $
     xz = qromb('cosm_intxz', zmin, z, /double)

  return, xz   ; Units are unknown

end

