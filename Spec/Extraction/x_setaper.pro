;+ 
; NAME:
; x_setaper   
;       Version 1.1
;
; PURPOSE:
;    Determine an aperture for an object automatically
;
; CALLING SEQUENCE:
;   x_setaper, spec, center, [frac], RADIUS=, /SKYSUB
;
; INPUTS:
;   spec      - Spectrum (fits file ok)
;   center    - center of aperture (guess, actual center not essential)
;   frac      - Fraction of profile to extract [default: 0.05 off each side]
;
; RETURNS:
;   aperture  - Aperture: float(2) array around center
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  RADIUS=  -- Region of object profile to consider
;  /SKYSUB  -- Calculate a median sky at outer regions of slit and 
;              subtract off
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_setaper, 'spec.fits', 200.
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Apr-2001 Written by JXP
;   02-May-2002 Revised sky level
;-
;------------------------------------------------------------------------------

function x_setaper, spec_in, center, frac, RADIUS=radius, SKYSUB=skysub

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'aper = x_setaper(spec, center, [frac], RADIUS=, /SKYSUB) [V1.0])'
    return, -1
  endif 

;  Optional Keywords

  ; Fraction of flux to throw away
  if not keyword_set( FRAC ) then frac = 0.05
  if not keyword_set( RADIUS ) then radius = 20L
  
; Read in fits file if necessary
  if size(spec_in,/type) EQ 7 then spec = xmrdfits(spec_in, /silent) $
  else spec = spec_in

; Setup
  npix = n_elements(spec)
  aper = fltarr(2)

;;;;;;;;;;;;;
; Find aperture

  ; Find the peak and set a region around it
  imx = round(center)
  rgmin = (imx - radius) > 0
  rgmax = (imx + radius) < (npix-1)
  ncen = rgmax-rgmin+1

  ; Subtract off the sky
  if keyword_set( SKYSUB ) then medlvl = 0. else begin
      medlvl = ( median(spec[rgmin:rgmin+4]) + median(spec[rgmax-4:rgmax]))/2.
  endelse

  ; Calculate the spline
  spl_x = findgen(ncen) + float(rgmin)
  spl_y = spec[long(spl_x)] - medlvl
  splin = spl_init(spl_x, spl_y, /double)

  ; Recalculate at 1000 points
  nspln = 5000L
  sval = fltarr(nspln)
  cumul = fltarr(nspln)
  tot = 0.
  dx = float(rgmax-rgmin)/float(nspln)
  for i=0L,nspln-1 do begin
      xval = float(rgmin) + dx*float(i)
      sval[i] = spl_interp(spl_x, spl_y, splin, xval, /double)
      if i GT 0 then cumul[i] = cumul[i-1] + sval[i]*dx $
      else cumul[i]=sval[i]*dx
      tot = tot+sval[i]*dx
  endfor

  ; Find edges
  mn = min(abs(cumul-tot*frac), imn)
  aper[0] = center-(float(rgmin) + dx*float(imn))
  mx = min(abs(cumul-tot*(1.-frac)), imx)
  aper[1] = (float(rgmin) + dx*float(imx))-center

  return, aper
end
