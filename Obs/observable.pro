;+ 
; NAME:
; observable
;   Version 1.1
;
; PURPOSE:
;    Calculate the airmass for an obj with a given DEC at an
;  observatory with a given DEC and an offset in RA
;
; CALLING SEQUENCE:
;  observable, month, obs_dec, obj_ra, obj_dec, PAR=par, $
;               POS=pos, L1=l1, L2=l2, KECK=keck, VLT=vlt
;
; INPUTS:
;  month   -- 
;  obs_dec -- DEC of the observatory
;  obj_ra  -- RA of the object
;  obj_dec -- DEC of the object
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /KECK  -- Calculate for Keck
;  /VLT   -- Calculate for VLT
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;written by E. Gawiser, last modified 05jul03
;anyone can use this code, please notify gawiser@astro.yale.edu regarding
;any errors found!    
;WARNING:  POS=0 is equivalent to not setting it, so use e.g. 0.001  
;WARNING:  "middle of the night" is not local midnight but is much more useful 
;in Chile it's about 00:40 in winter, 01:40 in summer with daylight savings!
;dec in degrees, RA in hours
;month is from 1.00 to 12.99, use decimals for days.  e.g automnal equinox
;is roughly 9.67.  (No accounting for 28/30/31 days)
;NOTE:  To get Keck Nasmyth platform limits, use /keck.  
;/par will output the amplitude of atmospheric dispersion from L1 to L2
;along the parallactic angle (default is 3000,10000A) and if POS is set it 
;will output the component perpendicular to a slit at position angle POS

;EXAMPLE FOR IMAGING:  
;observable, 1.7, -30., 3.5, -27. 
;models observing the CDF-S (03 30 -27) on January 21 from CTIO 
;The program outputs the RA overhead, length of night, and airmass values.

;EXAMPLE FOR SPECTROSCOPY:  
;observable, 10.3, 20., 3.5, -27., /par, /keck, l1=3500, l2=6000, pos=0.01
;models observing the CDF-S (03 30 -27) on October 9 from Keck with 
;slit(s) at position angle = 0.01
;The program outputs the Nasmyth platform limits, RA overhead, 
;length of night, airmass, parallactic angle, total amplitude of ADR from 
;3500 to 6000A, and amplitude of that ADR perpendicular to POS=0.01 
;for this field.  
;you can call /vlt instead of /keck to use atmospheric conditions at 
;Paranal (about 50% more ADR) instead of Mauna Kea.  Calling neither uses 
;atmospheric conditions at La Silla which cause a bit more ADR than Paranal.  
;Note that parallactic angle is 
;defined as angle towards zenith i.e. from blue to red light of the object,
;with North=0, East=90, etc. This has been tested versus the skycalc code 
;and works fine anywhere on the sky (doesn't mean much once the object is 
;below the horizon however ;) ).  
;------------------------------------------------------------------------
;
; EXAMPLES:
;  print, airmass( 20., 45., [-3, -2, -1, 0])
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Written by E. Gawiser
;-
;------------------------------------------------------------------------------
PRO observable, month, obs_dec, obj_ra, obj_dec, PAR=par, $
                POS=pos, L1=l1, L2=l2, KECK=keck, VLT=vlt

;leave in degrees and hours respectively before calculating airmass
time = findgen(13)-6. ;0 is middle of night, 6 hours each side, hour increments
ra_midnight = 2.*(month - 9.67)
if (ra_midnight lt 0.) then ra_midnight = ra_midnight + 24.
object_vector = findgen(n_elements(time))+1.  ;don't want 0.
object_vector = object_vector/object_vector   ;now all = 1.
object_vector = object_vector*(ra_midnight - obj_ra)
; that's the hour angle that the object is at at midnight, turned into 
; a vector the same dimensions as time
hour_angle = object_vector + time 

temp = where(hour_angle gt 12.)
if temp[0] ge 0 then hour_angle(where(hour_angle gt 12.)) = $
	hour_angle(where(hour_angle gt 12.)) - 24. 
temp = where(hour_angle le -12.)
if temp[0] ge 0 then hour_angle(where(hour_angle le -12.)) = $
	hour_angle(where(hour_angle le -12.)) + 24. 
;let's keep hour_angle between -12. and 12. for clarity although it 
;shouldn't matter!

;at Keck, 12 degree leads to 11.4 hrs Dec 21 -> 11.2 Jan -> 10.8 Feb ->
;10.3 Mar -> 9.7 Apr -> 9.1 May -> 8.9 hrs Jun 21 
;i.e. for 18 degree avg=9.2, amp=1.2 at abs(latitude)=20
;since 12-18 deg time is 27min Dec and 31 min Jun
;
;at Palomar,  12 degree gives 12.1 hrs Dec 21 -> 11.8 Jan -> 11.0 Feb ->
;10.1 Mar -> 8.9 Apr -> 7.9 May -> 7.5 Jun 21 
;i.e. for 18 degree avg= 8.9, amp=2.2 at abs(latitude)=34
;since 12-18 deg time is 30min Dec, 29 min Mar/Sep and 39 min Jun 
;the roughly sinusoidal variation has a larger amplitude in summer but 
;we don't need too much accuracy here so we'll keep it simple
;
;all likely observatories lie in this range of abs(dec) so should be ok
;to use linear interpolation for avg, amp:
night_avg = 9.2 - 0.3*(abs(obs_dec)-20.)/14.
night_amp = 1.2 + 1.0*(abs(obs_dec)-20.)/14. 

;angle to use is 3.14*(month-0.7)/6. which =0. at Dec 21, pi/2 at 
;Mar 21, pi at Jun 21, 3pi/2 at Sep 21, and 2 pi at Dec 21 again.
night_angle = 3.14*(month-0.7)/6. 

if obs_dec ge 0. then night_length = night_avg + night_amp * cos(night_angle) $
else night_length = night_avg - night_amp * cos(night_angle)
;reverse min/max in Southern hemisphere 

if keyword_set(KECK) then begin 
 print, 'Both Kecks vignetted by bottom shutter at airmass > 3.2'
 print, 'Keck I (LRIS/HIRES):  Nasmyth platform limit at airmass > 1.82 rising'
 print, '                  (except at dec < -30 or dec > 75)'
 print, 'Keck II (ESI/DEIMOS): Nasmyth platform limit at airmass > 1.67 setting'
 print, '                  (except at dec < -30 or dec > 60)'
 ;print, 'Current Keck I shutter vignetting at airmass <= 1.04'
 ;print, 'i.e. elevation >= 75 with maximum vignetting of half at zenith'
 ;shutter vignetting supposedly fixed now that both sides of it move :)
endif 

print, ' ' 
print, ra_midnight, ' is roughly overhead at the middle of the night'
print, night_length, ' is rough length of night between 18 deg. twilights'
print, ' '

obs_dec = obs_dec * 3.14159/ 180.
obj_dec = obj_dec * 3.14159/ 180.
hour_angle = hour_angle * 15. * 3.14159/ 180.
; can be an array of many hour angles like (findgen(13) - 6.)
; true angle is arccos of dot product of unit vectors pointing at (0,obs_dec)
; and (hour_angle, obj_dec) and airmass is sec(angle) so we just need
; airmass = 1/dot product 
; unit vectors given by (cos RA cos dec, sin RA cos dec, sin dec)
dummy = findgen(n_elements(hour_angle)) + 1.
dummy = dummy/dummy  ;should be all 1's
obs_vector = findgen(3,n_elements(hour_angle))
obs_vector(0,*) = cos(obs_dec)
obs_vector(1,*) = 0. 
obs_vector(2,*) = sin(obs_dec)
obj_vector = findgen(3,n_elements(hour_angle))
obj_vector(0,*) = cos(hour_angle)*cos(obj_dec)
obj_vector(1,*) = sin(hour_angle)*cos(obj_dec)
obj_vector(2,*) = dummy*sin(obj_dec)
cos_zenith_angle = total(obs_vector * obj_vector,1)

;airmass = 1. / cos(angle) ;= sec(angle)
airmass = 1/cos_zenith_angle

if keyword_set(PAR) then begin   ;this goes almost to the end of the code!

 sin_zenith_angle = sqrt(1. - cos_zenith_angle^2)

  ;use formula (9) from Filippenko 1982 based on law of sines on spherical 
  ;triangle from zenith to celestial pole to object:  
  ;sin parallactic / sin (90. - obs_dec) = sin (hour_angle) / sin(zenith_angle)

 sin_parallactic = sin(hour_angle) * sin(1.5708 - obs_dec) / sin_zenith_angle

 parallactic = asin(sin_parallactic) * 180./3.14159
	;want that in degrees for output, will be -90. to 90. 
	;negative hour angle gives negative parallactic, as desired 
	;since the spherical angle is from N to E for HA>0 but N to W for HA<0

 cos_parallactic = (cos_zenith_angle*sin(obj_dec)*(-1.) + sin(obs_dec))  / $
	sin_zenith_angle*cos(obj_dec)
  ;use law of cosines on sphere to figure out if abs(parallactic)>90. 

 temp = where(cos_parallactic lt 0.)
 if temp(0) ge 0 then $  ;i.e. this case applies
	parallactic(where(cos_parallactic lt 0.)) = $ 
	parallactic(where(cos_parallactic lt 0.))*(-1.) + 180.
  ;turns obtuse angles into 180. - acute_arcsin

 temp = where(parallactic gt 180.)
 if temp(0) ge 0 then $  ;i.e. this case applies
 	parallactic(where(parallactic gt 180.)) = $
	parallactic(where(parallactic gt 180.)) - 360. 
  ;puts it all into -180. to 180. range
  ;parallactic(where(parallactic lt -180.)) = $
  ;	parallactic (where(parallactic lt -180.)) + 360. 
  ;this one should be unnecessary, it's only the pi - P bit that makes 
  ;negative P into >180. 

;------------------------
;below code adapted from diff_atm_refr.pro, 
;     written by Enrico Marchetti, ESO, January 2001, 
;        found on ESO's website discussing VIMOS/La Silla atm.dispersion
; The original code computed a table of Differential Atmospheric Dispersion
; for a given zenithal distance "Z" for different wavelengths "Lambda"
; with respect to a reference wavelength "Lambda0".
;
; The atmospheric parameters can be adjusted to those characterstic
; of the site the computation is made for.
; The parameters listed below refer to the average La Silla conditions.
; (which should be appropriate for Pachon and Campanas, but not necessarily
; Mauna Kea - although the difference in output is so small between La Silla
; and Paranal that it's probably not a huge issue. - E. Gawiser July 4 2003)  
;
;--------------------------------------------------------
; wavelength range and reference wavelength to edit
;--------------------------------------------------------

 if keyword_set(l1) then lambda1=l1/10000. else begin 
    l1=3000.            ; [A]
    lambda1=l1/10000.   ; [um]
 endelse
 if keyword_set(l2) then lambda2=l2/10000. else begin
    l2=10000.           ; [A] 
    lambda2=l2/10000.   ; [um]
 endelse

;---------------------------------------------------------
; atmospheric parameters to edit
;---------------------------------------------------------
;default set to La Silla conditions (close enough for Pachon/LCO since 
;very similar results from Paranal values) - E. Gawiser 4jul03
;typical results compare well to those in Filippenko 1982
;results for Mauna Kea look about 30% low compared to NIR values found 
;by Gemini GNIRS team, not sure why since we adopted their TC and P values...
;also look 15% low compared to Drew Phillips' results plotted on proposed 
;LRIS ADC webpage - these cover full optical so are a more relevant check

if keyword_set(vlt) then begin  ;typical Paranal conditions (quite similar)
  print, 'Using typical atmospheric conditions for Paranal'
  TC=11.5D			;Temperature [C]
  RH=14.5D			;Relative Humidity [%]
  P=743.0D			;Pressure [mbar]
endif else begin
  if keyword_set(keck) then begin ;higher and colder at Mauna Kea!  
    print, 'Using typical atmospheric conditions for Mauna Kea'
    TC=0.0D			;Temperature [C]
    RH=50.0D			;Relative Humidity [%]
     ;currently a guess but it's nowhere near as dry as Paranal... 
     ;and varying this from 0.0 to 100.0 makes <1% change anyways 
    P=452.0D			;Pressure [mbar]
  endif else begin   ;the default
    print, 'Using typical atmospheric conditions for La Silla/Pachon/LCO'  
    TC=11.5D			;Temperature [C]
    RH=44.0D			;Relative Humidity [%]
    P=772.0D			;Pressure [mbar]
  endelse
endelse

;----------------------------------------------------------

 ZD = acos(cos_zenith_angle)   ;should be in radians
 T=TC+273.16D  ;now in Kelvin

 PS=-10474.0+116.43*T-0.43284*T^2+0.00053840*T^3

 P2=RH/100.0*PS
 P1=P-P2

 D1=P1/T*(1.0+P1*(57.90D*1.0E-8-(9.3250D*1.0E-4/T)+(0.25844D/T^2)))
 D2=P2/T*(1.0+P2*(1.0+3.7E-4*P2)*(-2.37321E-3+(2.23366/T)-(710.792/T^2)+(7.75141E4/T^3)))

 S0 = 1.0/Lambda2  ;used to be reference wavelength Lambda0
 S = 1.0/Lambda1   

 N0_1=1.0E-8*((2371.34+683939.7/(130-S0^2)+4547.3/(38.9-S0^2))*D1+$
     (6487.31+58.058*S0^2-0.71150*S0^4+0.08851*S0^6)*D2)

 N_1=1.0E-8*((2371.34+683939.7/(130-S^2)+4547.3/(38.9-S^2))*D1+$
     (6487.31+58.058*S^2-0.71150*S^4+0.08851*S^6)*D2)

 DR=Tan(ZD)*(N0_1-N_1)*206264.8

;end of "borrowed" code...  
;-----------------------

 amplitude = abs(DR)   
;no point showing a minus sign here since parallactic already identifies 
;the angle from blue to red.  That minus sign is flipped by swapping 
;l1 and l2 anyway.  

 if keyword_set(POS) then begin 
  angular_difference = parallactic - POS
  DR_perp = DR * sin(angular_difference * 3.14159/180.)
  amplitude_perp = DR_perp
;we want sign information here to remind us that blue and red are 
;swapping sides of the slit as meridian is crossed - this is important
;especially in comparing a science wavelength with the guider wavelength
;that difference is effectively set to zero during slitmask alignment 
;(for whatever effective wavelengths the imaging filter and guider filter have)
;but then as the field rotates only the guider light is kept constant!!!  
;we could use abs(DR_perp) and let the sign of parallactic reflect meridian
;crossing but for now this seems preferable

  print, 'For hours versus middle of night here is '
  print, '             airmass, parallactic angle (blue to red), amplitude (")'
  print, ' of differential refraction from', l1, ' to ', l2, ' Angstroms'   
  print, ' 1st amplitude is total, 2nd is perpendicular to POS=', POS
  for i=0,n_elements(airmass)-1 do print, time(i), airmass(i), parallactic(i),$
       amplitude(i), amplitude_perp(i) 
 endif else begin ;case of no POS specified
  print, 'For hours versus middle of night here is'
  print, '             airmass, parallactic angle (blue to red), amplitude (")'
  print, ' of differential refraction from', l1, ' to', l2, ' Angstroms'   

  for i=0,n_elements(airmass)-1 do print, time(i), airmass(i), parallactic(i),$
       amplitude(i) 
  ;
 endelse 
endif else begin  ;case of no /par requested 
  	;
  	print, 'For hours versus middle of night here is airmass:'
  	for i=0,n_elements(airmass)-1 do print, time(i), airmass(i)
endelse


;stop
end


