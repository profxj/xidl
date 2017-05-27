;+
; NAME:
;   discrete_arcfit_guess
;
; PURPOSE:
;   Find initial wavelength guess using discrete cross-correlation 
;
; CALLING SEQUENCE:
;   tset = discrete_arcfit_guess(spec, lambda, color=color, grating=grating, $
;                 func=func)
;
; INPUTS:
;   spec    - the spectrum 
;   lambda  - line list to match to
;   grating - grating (1200 or 831) [lines/mm]
;   lambda_c     -- central wavelength [Ang]
;   color   - 'red' or 'blue' side of DEIMOS
; 
; KEYWORDS:
;   func    - functional form to use
;
; OUTPUTS:
;   tset    - traceset mapping lines found in spec onto lambda list. 
; 
; COMMENTS:
;   This is much faster than lin_arcfit_guess and does not depend on 
;     any external optical models
;   --->  Should use real saturation/bad pixel mask
;
; REVISION HISTORY:
;
;       Thu Apr 25 17:12:33 2002, Doug Finkbeiner (dfink)
;		Written 
;
;----------------------------------------------------------------------
function discrete_arcfit_guess, spec, lambda, grating, lambda_c, $
                color=color, func=func, plot=plot
;, arcsat=arcsat

  if NOT keyword_set(color) then message, 'must set color!'
  if NOT keyword_set(grating) then message, 'grating keyword must be set!'
  if NOT keyword_set(lambda_c) then message, 'central wavelength not set!'
  
; to deal with refined 831/1200 lines/mm numbers
  grating=floor(grating)

; there may be multiple pixels flagged on top of saturated peaks!
; these are useful for initial guess; not for final fit. 
  
;  if n_elements(arcsat) eq 0 then sat   = (spec GT 60000) else $
;        sat=arcsat   ; replace this with a real saturation mask


; since we're doing discrete_correlate, we probably want to keep
 ; saturated peaks ????  JAN 1/15/02
 ;  peaky = findpeaks(spec, (1B-sat), nsig=10)
   peaky=findpeaks(spec,nsig=6)
 
  xpeak = where(peaky)

  if n_elements(xpeak) LE 1 then begin 
     message, 'not enough peaks!', /info
     print, 'not enough peaks!'
     return, 0
  endif else print,n_elements(xpeak),' peaks found'

  lagrange=[-500, 500]*1200./grating  ; in Angstroms (same units as lambda)
  if grating eq 1200 then begin 
     if color eq 'blue' then acoeff = [lambda_c-661., 661., 2., 0.,0.,0.]
     if color eq 'red'  then acoeff = [lambda_c+643., 643., -6, 0.,0.,0.]
  endif 
  if floor(grating) eq 831 then begin 
     ; don't know about these...
     pixscale = .460 * 2048.

     acoeff = [lambda_c-pixscale, pixscale, 5.74377, -0.83365,0.,0.] ; 831? l/mm grating
     if color eq 'red' then acoeff = [lambda_c+pixscale, pixscale, -5.4546, -0.9348,0.,0.] 
  endif 
  if grating eq 900 then begin
;     if color eq 'blue' then acoeff = [lambda_c-881, 0.40*2048, 5.,
;     0.]

; for Tomasso:

     if color eq 'blue' then acoeff = [lambda_c-865.3, 865.3, 7.8, 0.,0.,0.]
     if color eq 'red' then acoeff = [lambda_c+880, 880., -4., 0.,0.,0.]
;7.E-7*(2048)^2, -3.E-10*(2048)^3]
  endif

  if grating eq 600 then begin
;	if color eq 'blue' then acoeff = [lambda_c-1312.,1312., -8., 10.]
	if color eq 'blue' then acoeff = [lambda_c-1290.,1290., 13.8, -1.6,0.,0.]
        if color eq 'red'  then acoeff = [lambda_c+1320.,1320., -3.5, -1.4,0.,0.]
  endif

  if NOT keyword_set(acoeff) then begin 
     print, 'Grating: ', grating
     message, 'What grating are you using?!?'
  endif 

  dcoeff = [0, 4096.*.01/2, 0., 0.,0.,0.] ;*1200./grating ; cover this range +/-
  
  nstep = round(21*1200./grating)
  xrange = [0, n_elements(spec)]

  aset = discrete_tset_match(xpeak, lambda, xrange, func=func, $
             lagrange=lagrange, plot=plot, $
             acoeff=acoeff, dcoeff=dcoeff, nstep=nstep, box=5*1200./grating)
  
  return, aset
end 









