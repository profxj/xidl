;+
; NAME:
;   long_surf_trace2d
;
; PURPOSE:
;  Fit a 2D surface to the output from trace_crude.  Useful for
;  tracing arc lines among other things.
;
; CALLING SEQUENCE:
;  long_surf_trace2d, t, y, xcen, xerr, surffit
;
; INPUTS:
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS: 
;
; COMMENTS:
; 
; EXAMPLES:
;
; BUGS:
;    
; PROCEDURES CALLED:
; 
;
; REVISION HISTORY:
;   27-May-2005 Written by J. Hennawi (UCB)
;-
;------------------------------------------------------------------------------
pro long_surf_trace2d, t, y, xcen, xerr, surffit, $
  nycoeff=nycoeff, ntcoeff=ntcoeff, $
  res=res, sigrej=sigrej, mask=mask, adderr=adderr

   if NOT keyword_set(sigrej) then sigrej=4.
   if N_ELEMENTS(adderr) EQ 0 then adderr=0.03

   ; xcen and xerr come from trace_crude

   npix   = n_elements(t)
   if NOT keyword_set(nycoeff) then nycoeff = 6
   if NOT keyword_set(ntcoeff) then ntcoeff = 4

   invvar = 1.0/(abs(xerr)+adderr)^2 * (xerr LT 90.0 AND xerr GT 0)

   work2d = dblarr(npix,nycoeff*ntcoeff)
   worky = flegendre(y[*], nycoeff)
   workt = flegendre(t[*], ntcoeff)

   for i=0,ntcoeff-1 do begin
     for j=0,nycoeff-1 do begin
        work2d[*,j*ntcoeff+i] = worky[*, j] * workt[*,i]
     endfor
   endfor

   mask = invvar GT 0
   for iiter=1,3 do begin
     work2disig = work2d * $
            ((sqrt(invvar*mask))[*] # replicate(1,ntcoeff*nycoeff))
     alpha = transpose(work2disig) # work2disig
     beta = transpose(work2disig) # ((xcen * sqrt(invvar*mask))[*])

     svdc, alpha, w, u, v, /double
     res = svsol(u, w, v, beta, /double)

     ;  choldc, alpha, p, /double
     ;  res = cholsol(alpha,p,beta)

     surffit = xcen * 0.0
     surffit[*] = work2d # res
 
     residual = (xcen-surffit)*sqrt(mask*invvar)
     outlier = where(residual^2 GT sigrej^2, noutlier)
     if noutlier EQ 0 then break

     mask[outlier] = 0
   endfor

return
end

