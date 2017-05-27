;+
; NAME:
;   discrete_correlate_match
;
; PURPOSE:
;   Find match between discrete sets with discrete_correlation
;
; CALLING SEQUENCE:
;   discrete_correlate_iter, x1, x2, step=step
;
; INPUTS:
;   x1   - array of numbers we wish to correlate with x2
;   x2   - array, possibly large; perhaps only a small subset match
;           with x1. 
;
; KEYWORDS:
;   step - step size for search
;   lag  - lagrange for search(passed to discrete_correlate)
;
; RETURNS:
;   ind - array of indices for x2, i.e. x1 matches x2[ind]
;
; EXAMPLES:
;
; COMMENTS:
;   We should make this handle the 2-D case (see Hogg's
;   offset_from_pairs routine). 
;   
; REVISION HISTORY:
;
;       Thu Feb 21 22:19:49 2002, Douglas Finkbeiner (dfink)
;		Written
;
;----------------------------------------------------------------------
function discrete_correlate_match, x1, x2_in, step=step, lag=lag

  n1 = n_elements(x1) 

; -------- PASS 1: get offset between x1 and x2
  nbest = long(n_elements(x1) * .85)
  x2 = x2_in - discrete_correlate(x1, x2_in, step=step, nbest=nbest,lag=lag)

; for each x1, get index of nearest x2
  ind = lonarr(n1)
  for i=0, n1-1 do ind[i] = where(min(abs(x1[i]-x2)) EQ abs(x1[i]-x2))

  dx = x1-x2[ind]


; -------- PASS 2: remove linear trend (i.e. adjust scale)
;  poly_iter, x1, dx, 1, 3, xfit, coeff=coeff

    onexy = transpose([[fltarr(n_elements(x1))+1], [x1]])
     hogg_iter_linfit, onexy, dx, dx*0.+1., coeff

;	coeff=[reform(coeff),0.]


  scale = 1+coeff[1]
  x2 = x2*scale

  xoffs = discrete_correlate(x1, x2, step=step, nbest=nbest,lag=lag)
  x2 = x2-xoffs

  for i=0, n1-1 do ind[i] = where(min(abs(x1[i]-x2)) EQ abs(x1[i]-x2))

  dx = x1-x2[ind]

; -------- PASS 3: tweak offset 
  print, 'Tweak offset', median(dx), ' pixels'
  x2 = x2+median(dx)
  for i=0, n1-1 do ind[i] = where(min(abs(x1[i]-x2)) EQ abs(x1[i]-x2))

  return, ind
end
