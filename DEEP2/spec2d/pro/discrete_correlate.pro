;+
; NAME:
;   discrete_correlate
;
; PURPOSE:
;   Find optimal lag between discrete sets by cross-correlation
;
; CALLING SEQUENCE:
;   discrete_correlate, x1, x2, step=step, nbest=nbest
;
; INPUTS:
;   x1   - array of numbers we wish to correlate with x2
;   x2   - array, possibly large; perhaps only a small subset match
;           with x1. 
;
; KEYWORDS:
;   step     - step size for search
;   lagrange - lag range (Not LaGrange!) to explore in same units as x1,x2
;   box      - box size to call a "match"
;   nbest    - use only "nbest" best matches for statistics. 
;
; KEYWORD OUTPUTS:
;   sdev     - sum of squares of relative offsets of the lists
;   nmatch   - number of elements of x1 that "matched" (within "box")
;               of an element of x2
;   ind      - indices (in x2) for matches with elements of x1. 
;                Value of -1 means no match was found. 
; RETURNS:
;   best lag (i.e. x1+lag ~= x2)
;
; EXAMPLES:
;
; COMMENTS:
;   We should make this handle the 2-D case (see Hogg's
;   offset_from_pairs routine). 
;   
; REVISION HISTORY:
;
;       Wed Apr 24 22:01:27 2002, Doug Finkbeiner (dfink)
;		Massive improvements, added many keywords
;
;       Thu Feb 21 22:19:49 2002, Douglas Finkbeiner (dfink)
;		Written
;
;----------------------------------------------------------------------
function discrete_correlate, x1, x2, step=step, lagrange=lagrange, $
                 box=box, sdev=sdev, nmatch=nmatch, ind=ind, nbest=nbest

  if NOT keyword_set(step) then step = 1

  n1  = n_elements(x1)

; -------- Generate set of lags to probe

  if keyword_set(lagrange) then begin 
     xlag = findgen((lagrange[1]-lagrange[0])/step+1)*step+lagrange[0]
     min1 = min(x1, max=max1)
     wkeep = where((x2 GT min1+lagrange[0]) AND (x2 LT max1+lagrange[1]), ct)
     if ct LT 5 then begin 
        print, 'Working between ', min1+lagrange[0], ' and', max1+lagrange[1]
        print, 'Not enough lines to run!!!'
        sdev=1.E10
        ind=lonarr(n1)-1
        return,0.
     endif 
     x2trim = x2[wkeep]
  endif else begin 
     minx = min(x2, max=maxx)
     maxf = max(x1)
     xlag = findgen((maxx-minx+maxf)/step)*step+minx-maxf
     x2trim = x2
  endelse 
  x2trim = reform(x2trim, n_elements(x2trim))

; -------- Store results in sdev
  sdev = fltarr(n_elements(xlag))

; -------- Loop over lags
  for j=0, n_elements(xlag)-1 do begin 

;    for each element of x1, get closest value of x2
     x1lag = reform(x1+xlag[j], n1)
     join = [x1lag, x2trim]
     sind = sort(join)
     nj = n_elements(sind) 
     w1 = where(sind LT n1)
     offs = (join[sind[(w1+1)<(nj-1)]]-x1lag) < $
       (x1lag-join[sind[(w1-1)>0]])
     
     if keyword_set(box) then begin  
        offs = (offs < box) > (-box) ; limits on penalty for bad match
     endif 

     if keyword_set(nbest) then begin ; use only nbest best matches
        soffs = sort(abs(offs))
        offs = offs[soffs[0:(nbest-1) < (n1-1)]]
     endif
     sdev[j] = total(offs^2) ; big if match is bad
  endfor

  best = mean(where(sdev EQ min(sdev))) ; average ind

  sdev = sdev[best]

; -------- return match indices 

  if arg_present(ind) or arg_present(nmatch) then begin 
     ind = lonarr(n1)
     for i=0, n1-1 do begin
        junk = min(abs((x1[i]+xlag[best])-x2), indi)
        ind[i] = indi
     endfor
     if keyword_set(box) then begin  
        offs = x1-x2[ind]+xlag[best]
        bad = where(abs(offs) GT box, nbad)
        if nbad gt 0 then ind[bad] = -1
        nmatch = n_elements(offs) - nbad
     endif 

  endif 

  return, xlag[best]
end
