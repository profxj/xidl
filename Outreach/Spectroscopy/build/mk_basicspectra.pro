@outr_basicspectra
resolve_all
;; 
print, 'Try all modes after this launches..'
wait, 2
outr_basicspectra

save, /routines, file=getenv('XIDL_DIR')+'/OUTR_BASICSPECTRA.SAV'

cd, getenv('XIDL_DIR')
spawn, 'tar -zcvf basicspectra.tar.gz OUTR_BASICSPECTRA.SAV Outreach/Spectroscopy/Images/quasar.jpg Outreach/Spectroscopy/Images/sun.jpg Outreach/Spectroscopy/Images/sodium_lamp.jpg Outreach/Spectroscopy/Images/sdss_spiral.jpg Outreach/Spectroscopy/Spectra/sdss_spiral.fit.gz Outreach/Spectroscopy/Spectra/vanden.fits'
spawn, '\mv basicspectra.tar.gz Outreach/Spectroscopy/build/'
spawn, '\rm OUTR_BASICSPECTRA.SAV'
cd, getenv('XIDL_DIR')+'/Outreach/Spectroscopy'


end
