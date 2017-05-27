;+
; NAME:
;   discrete_correlate2d
;
; PURPOSE:
;   Find optimal lag between discrete sets by cross-correlation, in 2-dimen.
;
; CALLING SEQUENCE:
;   discrete_correlate, x1, x2, step=step, nbest=nbest
;
; INPUTS:
;   x1   - two-vector array of numbers we wish to correlate with x2
;   x2   - two-vector array, possibly large; perhaps only a small subset match
;           with x1. 
;
; KEYWORDS:
;   step     - step size for search, same in each dimension
;   lagrange - lag range (Not LaGrange!) to explore in same units as x1,x2
;   box      - box size to call a "match"
;   nbest    - use only "nbest" best matches for statistics. 
;   tmax     - maximum angular rotation explored (+-)
;   tstep    - size of angular step in search
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
;       Thu Jun 13,  Marc Davis, added 2-d capability
;
;       Wed Apr 24 22:01:27 2002, Doug Finkbeiner (dfink)
;		Massive improvements, added many keywords
;
;       Thu Feb 21 22:19:49 2002, Douglas Finkbeiner (dfink)
;		Written
;
;----------------------------------------------------------------------
function discrete_correlate2d, x1, x2, step=step, lagrange=lagrange, $
                 box=box, sdev=sdev, nmatch=nmatch, ind=ind, nbest=nbest, $
                  tmax=tmax, tstep=tstep

  if NOT keyword_set(step) then step = 10
  if NOT keyword_set(lagrange) then lagrange = [-10., 10.] ;10" lag is default

  
  n1  = (size(x1, /dimen))[0]

; -------- Generate set of lags to probe

;  if keyword_set(lagrange) then begin 
     xlag = findgen((lagrange[1]-lagrange[0])/step+1)*step+lagrange[0]
     ylag = xlag ;lags for 2nd dimension
     tlag = (tstep gt 0) ? findgen(2.*tmax/tstep +1)*tstep -tmax : 0

     x2trim = x2 ;keep all
; -------- Store results in sdev
  sdev = fltarr(n_elements(xlag), n_elements(ylag), n_elements(tlag))

; -------- Loop over lags

 
  x1lag = x1*0.
  x1r = x1*0.
  nxlags = n_elements(xlag)
  nylags = n_elements(ylag)
  for jx=0, nxlags-1 do begin ;loop over x
  for jy=0, nylags-1 do begin ;loop over y
  for jtheta=0, n_elements(tlag)-1 do begin ;loop over theta

     theta = tlag(jtheta)
     ctheta = cos(theta)
     stheta = sin(theta)
     x1r[*, 0]=  x1[*, 0]*ctheta + x1[*, 1]*stheta
     x1r[*, 1]= -x1[*, 0]*stheta + x1[*, 1]*ctheta

     offs =  x1r[*, 0]*0.
     x1lag[*, 0] = x1r[*, 0]+xlag[jx]
     x1lag[*, 1] = x1r[*, 1]+ylag[jy]

     for i=0,n1-1 do $
       offs[i]=min((x2[*, 0]- x1lag[i, 0])^2 + (x2[*, 1]- x1lag[i, 1])^2)

     if keyword_set(box) then begin  
        offs = (offs < box^2)  ; limits on penalty for bad match
     endif 
 
     if keyword_set(nbest) then begin ; use only nbest best matches
        soffs = sort(offs)
        offs = offs[soffs[0:(nbest-1) < (n1-1)]]
     endif
     sdev[jx, jy, jtheta] = total(offs) ; big if match is bad
  endfor
  endfor
  endfor

  ibest =where(sdev EQ min(sdev)) ; average ind


; -------- return match indices 

;  if arg_present(ind) or arg_present(nmatch) then begin 
;     ind = lonarr(n1)
;     for i=0, n1-1 do begin
;        junk = min(abs((x1[i]+xlag[best])-x2), indi)
;        ind[i] = indi
;     endfor
;     if keyword_set(box) then begin  
;        offs = x1-x2[ind]+xlag[best]
;        bad = where(abs(offs) GT box, nbad)
;        if nbad gt 0 then ind[bad] = -1
;        nmatch = n_elements(offs) - nbad
;     endif 
;
;  endif
 
  return, [xlag[ibest mod nxlags], ylag[fix(ibest/nxlags) mod nylags], $
        tlag[fix(ibest/nxlags/nylags)]]
end














