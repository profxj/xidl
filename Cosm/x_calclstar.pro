;+ 
; NAME:
; x_calclstar
;  Version 1.1
;
; PURPOSE:
;  Calculate apparent magnitude for Lstar at a given 
;      redshift for a given cosmology and Hubbles constant.
;
; CALLING SEQUENCE:
;   lstar = x_calclstar( z, H0=, AMIN=, /SILENT, OM=, OV=)
;
; INPUTS:
;    z -- Redshift
;
; RETURNS:
;    lstar -- Apparent magnitude for Lstar
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   lstr0 -- Assumed value of Lstar at z=0  [default: -21.06]
;   H0 -- Hubbles constant (km/s/Mpc)
;   OM -- Omega Dark Matter
;   OV -- Lambda
;   /SILENT -- No printing to the screen
;
; OPTIONAL OUTPUTS:
;  AMIN -- Mpc per arcmin at z [physical!]
;
; COMMENTS:
;
; EXAMPLES:
;   lstar = x_calclstar(2., H0=75., AMIN=amin)
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   22-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

function x_calclstar, z, H0=h0, AMIN=AMIN, SILENT=silent, $
		OM=om, OV=ov, LSTR0=lstr0, _EXTRA=extra

  common cosmolgy_cmmn, cosm_dm, cosm_K, cosm_h, cosm_Ob, $
     cosm_L, cosm_r, sigma_8

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'lstar = x_calclstar(z, H0=, AMIN=, /SILENT, '+$
	     	'LSTR0=, OM=, OV=) [v1.1]'
    return, -1
  endif 

;  if not keyword_set( LSTR0 ) then LSTR0 = -20.4  ;; Old value
  if not keyword_set( LSTR0 ) then LSTR0 = -21.12  ;; Blanton et al. 2003 (R band)

  ;; Initialize the common block
  cosm_common, H0=h0, Omegavac=OV, OmegaDM=OM, SILENT=silent, _EXTRA=extra
  if not keyword_set(silent) then $
    print, 'x_calclstar: Assuming L_0 = ', strtrim(LSTR0,2), $
    ' and H_0 = ', strtrim(cosm_h,2)

  ;; distance
  r1 = cosm_dist(z, SILENT=silent)
  DA = (!pi/180./3600.)*r1/(1+z)         ; Mpc/arcsec (physical)
  dL = r1*(1+z)                          ; Mpc (luminosity distance)

  ;; Lstar (apparent)
  lstar = 5*alog10(dL*1e5) + LSTR0
  
  ;; LGTH
  amin = 60*DA  ; Mpc/arcmin
  
  ;; Print
  if not keyword_set(silent) then begin
	print, 'x_calclstar:  Mpc per arcmin (physical) = ', amin
	print, 'x_calclstar:  Lstar = ', lstar
  endif	


  return, lstar

end

