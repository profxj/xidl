FUNCTION airmass, obs_dec, obj_dec, hour_angle ;angles in degrees except for HA

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'am = airmass(obs_dec, obj_dec, [hour_angle]) [V1.0]'
    return, -1
  endif 

; Optional keywords

if not keyword_set( HOUR_ANGLE ) then hour_angle = findgen(5) - 2.




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
cos_angle = total(obs_vector * obj_vector,1)
;angle = acos(cos_angle)
;airmass = 1. / cos(angle) ;= sec(angle)
airmass = 1/cos_angle


;stop
return, airmass
;stop
end

;WARNING:  for reasons that I don't understand, the input variables
;get modified and returned, even if their names when sent differ 
;from the names used in this code.  Hence repeated calls will gradually 
;alter the array of hour_angle.  
