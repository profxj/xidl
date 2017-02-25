FUNCTION LONG_FLUXNORM, wave, flux, sig, mag, filter, MED_WIDTH = MED_WIDTH

IF NOT KEYWORD_SET(MED_WIDTH) THEN MED_WIDTH = 10L
; f_lambda= 1.0e-17 erg/cm^2/s/A = 1.0e-17*1.0e8 erg/cm^2/s/cm
; 1 Jansky = 1e-23 erg/cm^2/s/Hz
c = 2.9979246d10
angstrom = double(1.0d-8)
jans_scale =  9.1865627d-11 ; 1.0e-17*1.0e8*(1000.0*angstrom)^2/c/Jansky/3631
LN10 = alog(10.0D)

f_lam_med = djs_median(flux, width = MED_WIDTH, boundary = 'reflect')
f_nu = reverse((wave/1000.0D)^2*f_lam_med)
lognu = reverse(alog10(c/(wave*angstrom)))
F_m = sdss_filter_nu(lognu, filter)
mag_int = int_tabulated(lognu, LN10*f_nu*F_m)
A_norm = (10.d^(-0.4D*mag))/(jans_scale*mag_int)

RETURN, A_NORM
END
