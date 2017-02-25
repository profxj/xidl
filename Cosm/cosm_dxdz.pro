;+ 
; NAME:
; cosm_dxdz
;
; PURPOSE:
;    Calculate dX/dz at a given redshift
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;   z -- Redshift
;
; RETURNS:
;   dX/dz -- Cosmological pathlength
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
  

function cosm_dxdz, z, zmin=zmin, OM=om, OV=ov, NOINIT=noinit, $
                  exact=exact, _EXTRA=extra

  common cosmolgy_cmmn

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'xz = cosm_xz(z, zmin= OM=, OV=, /W05MAP, /NOINIT) [v1.1]'
    return, -1
  endif 
  
  if not keyword_set(zmin) then zmin = 0.

  if not keyword_set(NOINIT) then $
     cosm_common, H0=h0, Omegavac=OV, OmegaDM=OM, _EXTRA=extra

  dXdz = (1.+z)^2 / sqrt(cosm_L + (cosm_K)*(1+z)^2 + cosm_dm*(1+z)^3)

  return, dXdz   ; Unitless

end

