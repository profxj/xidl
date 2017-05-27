;+
; NAME:
;   ivarsmooth
;
; PURPOSE:
;   to generate inverse variance smoothed 1d spectrum
;
; CALLING SEQUENCE:
;   spec = ivarsmooth(flux,ivar,window,[outivar])
; 
; INPUTS:
;   flux -- spectrum to be smoothed
;   ivar -- inverse variance of this spectrum
;   window -- smoothing length in pixels (turned into odd number if
;             even or float)
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;    spec -- resulting smoothed 1d array
;
; OPTIONAL OUTPUTS:
;   outivar -- resulting inverse variance
;
; RESTRICTIONS:
;    
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;  JAN 04jul2002, comments by MD
;----------------------------------------------------------------------

function ivarsmooth, flux, ivar, window, outivar

  s = size(flux)
  ndim = s[0]
  nflux = n_elements(flux)
  simpleflux = ndim eq 1 OR ndim eq 2 and (s[2] > s[1]) eq nflux

  if NOT simpleflux then begin
     message, 'Flux must be a 1-D array!', /continue
     return, flux
  endif

  flux = reform(flux)
  ivar = reform(ivar)

  halfwindow = floor( (round(window)-1)/2)

  shiftarr = fltarr(nflux, 2*halfwindow+1)
  shiftivar = shiftarr
  shiftindex = shiftarr
  indexarr = findgen(nflux)
  indnorm = (fltarr(2*halfwindow+1)+1) ## indexarr

  for i=-halfwindow, halfwindow do begin & $
     shiftarr(*, i+halfwindow) = shift(flux, i) & $
     shiftivar(*, i+halfwindow) = shift(ivar, i) & $
     shiftindex(*, i+halfwindow) = shift(indexarr, i) & $
  endfor

  wh = where( abs(shiftindex-indnorm) gt (halfwindow+1))
  shiftivar[wh] = 0.

  outivar = total(shiftivar, 2)
  nzero =  where(outivar gt 0,zeroct)
  smoothflux = total(shiftarr*shiftivar, 2)
  if zeroct gt 0 then smoothflux[nzero] = smoothflux[nzero]/outivar[nzero] $
     else smoothflux = smooth(flux,2*halfwindow+1);kill off NAN's

return, smoothflux
end
