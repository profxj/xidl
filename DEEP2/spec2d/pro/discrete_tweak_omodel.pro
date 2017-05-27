;+
; NAME:
;   discrete_tweak_omodel
;
; PURPOSE:
;   Refine initial wavelength guess using discrete cross-correlation 
;
; CALLING SEQUENCE:
;   tset = discrete_tweak_omodel(spec, lambda, coeff,  $
;                 func=func)
;
; INPUTS:
;   spec    - the spectrum 
;   lambda  - line list to match to
;   coeff   - estimate of coefficients from optical model
; 
; KEYWORDS:
;   func    - functional form to use
;
; OUTPUTS:
;   tset    - traceset mapping lines found in spec onto lambda list. 
; 
; COMMENTS:
;   This is much faster than discrete_arcfit_guess, but does not depend on 
;     an external optical model
;   --->  Should use real saturation/bad pixel mask
;
; REVISION HISTORY:
;
;       Thu Apr 25 17:12:33 2002, Doug Finkbeiner (dfink)
;		Written 
;
;----------------------------------------------------------------------
function discrete_tweak_omodel, spec, lambda, coeff, func=func, plot=plot
;, arcsat=arcsat
  

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
       peaky=findpeaks(spec,nsig=4)
       xpeak = where(peaky)
  endif

  if n_elements(xpeak) LE 1 then begin 
     message, 'not enough peaks!', /info
     print, 'not enough peaks!'
     return, 0
  endif 
;else print,n_elements(xpeak),' peaks found'

  linscale=round(coeff[1]/650.*4)/4.  
; ratio of angstroms/pix compared to 1200, to nearest 0.25

  lagrange=[-20., 20.]*linscale  ; in Angstroms (same units as lambda)

  if n_elements(coeff) eq 0 then message, 'Optical model estimate required!'

  acoeff=coeff

  dcoeff = [0, 2.5, 0., 0.,0.,0.]*linscale ; cover this range +/-
  
  nstep = 15
  xrange = [0, n_elements(spec)]

  aset = discrete_tset_match(xpeak, lambda, xrange, func=func, $
             lagrange=lagrange, plot=plot, $
             acoeff=acoeff, dcoeff=dcoeff, nstep=nstep, $
             stepscale=0.1*linscale,box=5) ;*linscale)

;nstep=11

  aset = discrete_tset_match(xpeak, lambda, xrange, func=func, $
             lagrange=lagrange/8., plot=plot, $
             acoeff=aset.coeff, dcoeff=dcoeff/15., nstep=nstep, $
             stepscale=0.006*linscale,box=5) ;*linscale)



  
  return, aset
end 









