PRO observable, month, obs_dec, obj_ra, obj_dec ;dec in degrees, RA in hours
;month is from 1 to 12, use decimals for days.(No accounting for 28/30/31 days)
;e.g. autumnal equinox is roughly 9.67.  

;leave in degrees and hours respectively before calling airmass
time = findgen(13)-6. ;0 is middle of night, 6 hours each side, hour increments
ra_midnight = 2.*(month - 9.67)
if (ra_midnight lt 0.) then ra_midnight = ra_midnight + 24.
object_vector = findgen(n_elements(time))+1.  ;don't want 0.
object_vector = object_vector/object_vector   ;now all = 1.
object_vector = object_vector*(ra_midnight - obj_ra)
; that's the hour angle that the object is at at midnight, turned into 
; a vector the same dimensions as time
hour_angle = object_vector + time 

print, 'Both Kecks vignetted by bottom shutter at airmass > 3.2'
print, 'Keck I (LRIS/HIRES):  Nasmyth platform limit at airmass > 1.82 rising'
print, '                  (except at dec < -30 or dec > 75)'
print, 'Keck II (ESI): Nasmyth platform limit at airmass > 1.67 setting'
print, '                  (except at dec < -30 or dec > 60)'
print, 'Current Keck I shutter vignetting at airmass <= 1.04'
print, 'i.e. elevation >= 75 with maximum vignetting of half at zenith'
print, ' ' 
print, ra_midnight, 'is roughly overhead at the middle of the night'
print, 'Now for hour increments versus the middle of the night of:'
print, time
print, 'The airmass of this object will be:'
print, airmass(obs_dec, obj_dec, hour_angle)


;stop
return
end

