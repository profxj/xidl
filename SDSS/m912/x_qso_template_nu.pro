;+ 
; NAME:
; qso_template_nu
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
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
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2005  Written by Joe Hennawi UCB 
;-
;------------------------------------------------------------------------------
FUNCTION X_QSO_TEMPLATE_NU, lognu

common kcor_common_nu, eigen_lognu, eigen_flux, eigen_flux_spline


IF NOT KEYWORD_SET(eigen_flux_spline) THEN BEGIN

;   vanden Berk template with host galaxy componenet removed
    new_file = getenv('XIDL_DIR')+'/SDSS/m912/dr1QsoSpec.gfix.fits'
    str = xmrdfits(new_file,1,/sile)
    lam_van = str.wav
    eigen_flux_lam = str.flux1
    eigen_flux2 = str.flux2
    delvarx, str
;    rdfloat, new_file, lam_van, eigen_flux_lam, eigen_flux2, skipline = 1
    c = 2.9979246d10
    angstrom = 1.0d-8 
    nu = c/(lam_van*angstrom)
    eigen_lognu = reverse(alog10(nu))
;   units are arbirtray here
    eigen_flux = reverse(1.0e-7*lam_van^2*eigen_flux_lam)
;   compute spline
    eigen_flux_spline = spl_init(eigen_lognu, eigen_flux)
ENDIF

min_log = min(eigen_lognu, i_min)
max_log = max(eigen_lognu, i_max)

ind_spline = WHERE(lognu GT min_log AND lognu LT max_log)
ind_min = WHERE(lognu LE min_log)
ind_max = WHERE(lognu GE max_log)

F_qso = dblarr(n_elements(lognu))

IF ind_spline[0] NE -1 THEN $
  F_qso[ind_spline] = $
  SPL_INTERP(eigen_lognu, eigen_flux, eigen_flux_spline, lognu[ind_spline]) 
IF ind_min[0] NE -1 THEN  F_qso[ind_min] = 0.0 
IF ind_max[0] NE -1 THEN  F_qso[ind_max] = 0.0 
RETURN, F_qso

END
