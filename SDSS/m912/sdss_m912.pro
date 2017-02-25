;+
; NAME:
;   sdss_m912
;
; PURPOSE:
;   Compute the apparent magnitude one would measure at the 
;   Lyman limit in the observed frame at lambda = (1 + z)912
;
; CALLING SEQUENCE:
;   m_912 = m912(z,m,filter,[ logLv = , ALPHA_EUV = , OMEGA0 = , LAMBDA0 = $
;                          , W = , LIT_H = ])
; INPUTS:
;   z         -  Redshift of the quasar
;   m         -  Apparent magnitude
;   filter    -  filter name SDSS [u,g,r,i,z,B,B_J]
;
; OPTIONAL KEYWORDS:
;   logLv      - Log_10 of the specific luminosity at the Lyman edge (ergs/s/Hz)
;   ALPHA_EUV  - Spectral slope in the extreme UV. Default is 1.57   
;                (Tefler et al 2002). A power law v^{-ALPHA_EUV} 
;                is pasted onto the quasar spectrum at 1285 A. 
;   OMEGA0     - Matter density. Default is 0.27
;   LAMBDA0    - Dark energy density. Default is 0.73
;   W          - Equation of state parameter. Default is w = -1.0
;   LIT_H      - Hubble h parameter. Required if logLv is specified
;
; OUTPUTS:
;   m_912       - Apparent magnitude one would measure at the Lyman 
;                 limit in the observed frame at lambda = (1 + z)912
;
; COMMENTS:
;   The code will return an error if one tries to use a filter with 
;   lambda_min < (1 + z)*1216. The lyman alpha forest will reduce the 
;   observed apparent magnitude, and the mean flux blueward of 1216 
;   is not reproduced correctly in the  quasar template used here. 
; 
; EXAMPLE:
;     Compute m_912 for a quasar at z=2.43 which has g=20.5
;     m_912=m912(2.43,20.5,1,logLv=logLv)
;     print,m_912
;       22.199551
;     print,loglv
;       29.812608
;
; PROCEDURES CALLED:
;   dofa()
;   sdss_filter_nu()
;   qso_template_nu()
;   obs_int_nu()
;   
; REVISION HISTORY:
;   17-Sep-2005  Written by Joe Hennawi UCB 
;-
;------------------------------------------------------------------------------
; Compute integrand of observed frame k-correction integral

FUNCTION OBS_INT_NU, lognu, z, F_m

; evaluate quasar flux at rest frame wavelength
nu_rest = (10.0D^lognu)*(1.0D + z)
lognu_rest = alog10(nu_rest)
F_qso = qso_template_nu(lognu_rest)
LN10 = alog(10.0D)

OBS = LN10*F_qso*F_m
RETURN, OBS
END


FUNCTION sdss_m912, znow, m, filter, logLv = logLv $
               , ALPHA_EUV = ALPHA_EUV, OMEGA0 = OMEGA0, LAMBDA0 = LAMBDA0 $
               , w = w, LIT_H = LIT_H, IGNORE = IGNORE

IF n_elements(znow) NE n_elements(m) THEN $
  message, 'The number of elements in znow and m must be the same'

; Fundamental constants we will need
IF NOT KEYWORD_SET(LIT_H) THEN LIT_H = 0.72
IF NOT KEYWORD_SET(OMEGA0) THEN OMEGA0 = 0.27
IF NOT KEYWORD_SET(LAMBDA0) THEN LAMBDA0 = 0.73
IF NOT KEYWORD_SET(W) THEN W = 1.0
; If no EUV spectral slope is given, use the Tefler 2002 value
IF NOT KEYWORD_SET(ALPHA_EUV) THEN ALPHA_EUV = 1.57D

HORIZON = double(2.9979246d9)   ; HORIZON in pc/h
parsec = double(3.086d18)
logHORpc = 27.966217
angstrom = 1.0d-8 
c = 2.9979246d10
nu_912  = c/(911.76D*angstrom)
nu_1285 = c/(1285.0D*angstrom)
nu_1285_right =  c/(1282.5D*angstrom)
nu_1285_left  =  c/(1287.5D*angstrom)
lognu_912   = alog10(nu_912)
lognu_1285  = alog10(nu_1285)
lognu_1285l = alog10(nu_1285_left)
lognu_1285r = alog10(nu_1285_right)
AB_source = 3631.0D   
; The standard AB magnitude reference flux is 
; 3631 Jy = 3631*1.0e-23 erg/cm^2/s/Hz

; This common block is needed for the integrals against the template
common kcor_common_nu, eigen_lognu, eigen_flux, eigen_flux_spline

; effective beginning wavelengths of SDSS filters
; these are defined as the wavelength where the 
; filters are 10% of their peak value
filter_names = ['u', 'g', 'r', 'i', 'z', 'B', 'B_J']
filter_ind = lonarr(n_elements(filter))
FOR j = 0L, n_elements(filter)-1L DO $
  filter_ind[j] = WHERE(filter_names EQ filter[j])
temp = WHERE(filter_ind EQ -1, nfilt_bad)
IF nfilt_bad NE 0 THEN message, 'Unidentified filters'
;            u        g       r      i        z      B    B_J
wave_beg = [3125.0, 3880.0, 5480.0, 6790.0, 8090.0, 3772.0, 3706.3]
; Since the template is unreliable in the Lya forest, require
; that the filter be redder than observed frame Lya
bad_inds = WHERE((1.0D + znow)*1215.67 GE wave_beg[filter_ind], nbad)
IF nbad NE 0 AND NOT KEYWORD_SET(IGNORE) THEN BEGIN
    splog, 'ERROR: The filter begins at ', wave_beg[filter_ind] $
           , ' which must be redder than observed frame Lya'
    splog, '       There are ', nbad, ' bad redshifts'
    stop
END

nz = n_elements(znow)
m_912 = dblarr(nz)
anow = 1.0/(1.0+znow)
; Compute Luminosity distance in cm
;DL1 = (1.0+znow)*dofa(anow, omega0, lambda0, w)/lit_h

;DL1 = cosm_dist(znow,/init,h0=lit_h*100.);*parsec*10e6
DL1=(1.0+znow)*comdis(znow,omega0,lambda0)/lit_h
logDL = logHORpc + alog10(DL1)

; create a vector of frequencies in log units
n_nu = 25000
dlognu = 0.0001
lognubeg = 14.0
lognu = dindgen(n_nu)*dlognu + lognubeg
; read in spline of the filter which will be used onto these frequencies
F_m = sdss_filter_nu(lognu, filter)

; Tefler et al. 2002 found that the break between power laws 
; occurred between 1200 and 1300 A. The cleanset regino of the 
; quasar template in this range (i.e. far from emission lines) is 
; between SiII 1265.22 A and OI 1305.4 A. Thus we match the quasar
; EUV power law to the optical template at 1285 A, but we average
; over 5 A as the template is a bit noisy. 
avg_pix = WHERE(lognu GE lognu_1285l AND lognu LE lognu_1285r, navg)
F_van = qso_template_nu(lognu)
F_van1285 = total(F_van[avg_pix])/double(navg)
log_scale = alog10(F_van1285) + alpha_EUV*(lognu_1285-lognu_912)

; integrate the vanden Berk qso template against this filter
FOR j = 0L, nz-1L DO BEGIN
;   redshift the rest frame spectrum to the observed frame 
;   and integrate the redshifted spectrum over the filter curve
    obs = obs_int_nu(lognu, znow[j], F_m)
    obs_integ = int_tabulated(lognu, obs, /DOUBLE)
    m_912[j] = -2.5D*(log_scale - alog10(obs_integ)) + m[j]
ENDFOR
; Now convert m912 to a specific luminosity in Janskies
logLv = alog10(anow*4.0D*!dpi*AB_source) + 2.0D*logDL - 0.4D*m_912 -23.0D

RETURN, m_912
END
