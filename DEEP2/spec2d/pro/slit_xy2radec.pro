;+
;
; NAME
;      slit_xy2radec.pro
;
; PURPOSE
;      Compute the distance in RA and DEC between point A and point B
;      as measured on the detector. It work like this...you have two
;      objects on your 2-d slit and you know the pixel distance in the
;      spatial direction (pixel_distance) separating the two objects,
;      but you'd also like to know the distance separating them in RA
;      and DEC on the sky. This assumes that you know the position of
;      one or at the least you can make a good guess at the
;      declination of the two sources. 
;
; SYNTAX
;      slit_xy2radec, pixel_distance, declination, $
;                     pixel_scale=pixel_scale, slitpa=slitpa, $
;                     maskpa=maskpa, delta_ra=delta_ra, $
;                     delta_dec=delta_dec, /arcseconds, /arcminutes,
;                     /degrees, /radians
;
; INPUTS
;      pixel_distance = the distance between point A and point B on
;                       the detector in pixels in the spatial
;                       direction.
;      declination = the declination corresponding to point A.
;      pixel_scale = the pixel scale of the detector assumed to be in
;                    units of arcseconds per pixel.
;      slitpa = the position angle of the slit relative to North on
;               the sky with positive angles increasing towards the
;               East. This parameter in conjunction with the maskpa
;               parameter allows for the possibility of tilted slits
;               where the slit is tilted with respect to the spatial
;               direction on the detector. 
;      maskpa = the position angle of the spatial axis of the mask
;               relative to North on the sky with positive angle
;               increasing towards the East.
;      delta_ra = a variable that will be set equal to the distance
;                 between points A and B in RA on the sky.
;      delta_dec = a variable that will be set equal to the distance
;                  between points A and B in DEC on the sky.
;
; KEYWORDS
;      /arcseconds OR /arcminutes OR /degrees OR /radians
;      /arcseconds = if this keyword is set, then the output
;                    parameters delta_ra and delta_dec will be
;                    returned in units of arcseconds on the sky.
;      /arcminutess = if this keyword is set, then the output
;                     parameters delta_ra and delta_dec will be
;                     returned in units of arcminutes on the sky.
;      /degrees = if this keyword is set, then the output
;                 parameters delta_ra and delta_dec will be
;                 returned in units of degrees on the sky.
;      /radians = if this keyword is set, then the output
;                 parameters delta_ra and delta_dec will be
;                 returned in units of radians on the sky.
;
; OUTPUTS
;      delta_ra = the distance from point A to point B on the sky in
;                 right ascension. The units are specified by the user
;                 according to the /arcseconds, /arcminutes, /degrees,
;                 and /radians keywords. The default setting is to
;                 return delta_ra in units of arcseconds.
;      delta_dec = the distance from point A to point B on the sky in
;                  declination. The units are specified by the user
;                  according to the /arcseconds, /arcminutes, /degrees,
;                  and /radians keywords. The default setting is to
;                  return delta_ra in units of arcseconds.
;      Note that the distance on the sky between points A and B are
;      thus given by sky_dist = sqrt( delta_ra^2 + delta_dec^2 ).
;
; PROCEDURES CALLED
;
;
; EXAMPLES
;      Compute the distance in RA and DEC between point A and point B
;      which are separated on the detector by a pixel distance of 5
;      pixels in the spatial direction. The declination of point A is
;      known to be dec = 36.125 degrees. The pixel scale of the object
;      is 0.207 arcseconds per pixel. The distance in RA and DEC will
;      be given in arcseconds:
;        slit_xy2radec, 5.0, 36.125, pixscl=0.207, slitpa=5.0, $
;              maskpa=15.0, delta_ra=delra, delta_dec=deldec
;
; COMMENTS
;      None.
;
; HISTORY
;      Created February 17, 2003 by mcc.
;
;-

; -----------------------
; define the following two simple functions which allow the user to
; supply the position angle info in degrees rather than radians.
function cosd, x
  y = cos( x / !radeg )
  return, y
end

function sind, x
  y = sin( x / !radeg )
  return, y
end
; -----------------------

pro slit_xy2radec, pixel_distance, declination, $
                   pixel_scale=pixel_scale, $
                   slitpa=slitpa, maskpa=maskpa, $
                   delta_dec=delta_dec, delta_ra=delta_ra, $
                   radians=radians, arcseconds=arcseconds, $
                   arcminutes=arcminutes, degrees=degrees

; check that enough parameters were passed.
  if n_params() lt 2 then begin
      message, '(slit_xy2radec.pro) Incorrect number of parameters!'
  endif

; define more terse (eaiser to type) variable names.
  dec = declination
  pixdist = pixel_distance

; check if any of the keywords were set.
  if not(keyword_set(slitpa)) then slitpa = 0.0
  if not(keyword_set(maskpa)) then maskpa = 0.0
  if keyword_set(pixel_scale) then pixscl = pixel_scale $
  else pixscl = 1.0

  scale_factor = 1.0
  if keyword_set(arcseconds) then scale_factor = 1.0
  if keyword_set(arcminutes) then scale_factor = 60.0
  if keyword_set(degrees) then scale_factor = 3600.

; calculate the distance between A and B on the sky.
  if keyword_set(radians) then begin
      delta_dec = pixdist * pixscl * cos(maskpa) / cos(slitpa - maskpa)

      delta_ra = pixdist * pixscl / cos(slitpa - maskpa) * $
        sin(maskpa) / cos(dec)
  endif else begin
      delta_dec = pixdist * pixscl * cosd(maskpa) / $
        cosd(slitpa - maskpa) / scale_factor

      delta_ra = pixdist * pixscl / cosd(slitpa - maskpa) * $
        sind(maskpa) / cosd(dec) / scale_factor
  endelse

end

