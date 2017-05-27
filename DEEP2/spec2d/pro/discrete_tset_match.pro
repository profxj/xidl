;+
; NAME:
;   discrete_tset_match
;
; PURPOSE:
;   Find traceset which transforms one list of numbers to match another
;
; CALLING SEQUENCE:
;   tset = discrete_tset_match(x, y, xrange, func=func, lagrange=lagrange, $
;                 acoeff=acoeff, dcoeff=dcoeff, nstep=nstep, box=box)
;
; INPUTS:
;   x      - set of numbers to match to y
;   y      - another set of numbers
;   xrange - xrange of traceset to find (the polynomials will be
;              orthogonal on this domain)
;   func   - functions to use (default: legendre polynomials)
;   lagrange - range of lags to explore
;   acoeff - center of (polynomial coefficient) parameter space to explore
;   dcoeff - range to explore, plus and minus from acoeff to explore
;
; KEYWORDS:
;   nstep  - number of steps to take while exploring
;   box    - numbers closer than "box" are called a match
;   plot   - to make a plot
; OUTPUTS:
;    returns traceset transforming x to a subset of y
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;
;       Wed Apr 24 22:06:48 2002, Doug Finkbeiner (dfink)
;		Written
;
;----------------------------------------------------------------------
; find the traceset that maps x to elements of y (some may be missing)
function discrete_tset_match, x, y, xrange, func=func, lagrange=lagrange, $
                 acoeff=acoeff, dcoeff=dcoeff, nstep=nstep, box=box, $
                 plot=plot,stepscale=stepscale
  
  if NOT keyword_set(xrange) then message, 'must set xrange!'
  if NOT keyword_set(nstep)  then nstep = 21
  if NOT keyword_set(acoeff) then message, 'must set acoeff!'
  if NOT keyword_set(dcoeff) then message, 'must set dcoeff!'
  if NOT keyword_set(func)   then func = 'legendre'
  if NOT keyword_set(box)    then box = 5
  if n_elements(stepscale) eq 0 then stepscale=1.

  nacoeff = N_elements(acoeff)
  aset = tset_struc(func, nacoeff)
  aset.xmin = xrange[0]
  aset.xmax = xrange[1]
 
  sig  = fltarr(nstep)
  dy   = fltarr(nstep)
  step = findgen(nstep)/(nstep-1)*2-1

  for i=0, nstep-1 do begin 
     
     aset.coeff = acoeff+step[i]*dcoeff
     traceset2xy, aset, x, ty
     
     dy[i] = discrete_correlate(ty, y, step=5*stepscale, lagrange=lagrange, $
                                box=box, sdev=sdev)
     sig[i] = sdev
  endfor 
  foo = min(sig, imin)

  if sdev eq 1.E10 then return,aset

  aset.coeff = acoeff+step[imin]*dcoeff
  aset.coeff[0] = aset.coeff[0]+dy[imin]

  traceset2xy, aset, x, ty
  dy = discrete_correlate(ty, y, step=1*stepscale, lagrange=lagrange/5, $
                          box=box, nmatch=nmatch, ind=ind, sdev=sdev)
   
  aset.coeff[0] = aset.coeff[0]+dy

  w = where(ind NE -1)
  if keyword_set(plot) then begin 
     wset, 0
     plot, ty[w], y[ind[w]]-ty[w], ps=7, yr=[-5, 5], xtit='wavelength [Ang]'
  endif 
  return, aset
end
