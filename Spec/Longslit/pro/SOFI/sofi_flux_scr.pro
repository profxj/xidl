spawn, 'mkdir Combine'
path = '/Users/joe/DATA/SOFI_DATA/'
tellfiles = path + ['sci-SOFI_0072.fits', 'sci-SOFI_0076.fits']
type = 'B9'
V = 9.0D
sensfunc = 'Combine/HD137005_0_sens.fits'
niri_sensfunc, tellfiles, type, sensfunc, V = V, magfunc = magfunc $
               , loglam = loglam, flux = flux, ivar = ivar, /CHECK

scifiles = path + ['sci-SOFI_0123.fits', 'sci-SOFI_0123.fits', 'sci-SOFI_0131.fits', 'sci-SOFI_0131.fits']

objid = [1, 2, 1, 2]
;;
outfile = 'Combine/J2003-3251.fits'
niri_fluxcal, scifiles, sensfunc, loglam = loglam, flux = flux $
              , ivar = ivar, mask = mask, outfile = outfile, OBJID = objid $
              , ARRFLUX = FLUX_ARR, ARRIVAR = IVAR_ARR $
              , SIGFLUX = flux_sig, SIGIVAR = ivar_sig $
              , /CHECK


END
