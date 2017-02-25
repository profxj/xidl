;+ 
; NAME:
; x_trace2d   
;     Version 1.2
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; 
; OUTPUTS TO SNGL:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   ??-2004 Written by SB
;                                
;------------------------------------------------------------------------------
pro x_trace2d, xcen, xerr, surffit, nycoeff=nycoeff, ntcoeff=ntcoeff, $
            res=res, sigrej=sigrej, t=t, y=y, mask=mask

  if  N_params() LT 3  then begin 
      print,'Syntax:  ' + $
        'x_trace2d(xcen, xerr, surffit, SIGREJ=, NTCOEFF=, NYCOEFF= [v1.2]'
      return
  endif 
  
   if NOT keyword_set(sigrej) then sigrej=4.

   ; xcen and xerr come from trace_crude

   ntrace = (size(xcen))[2]
   npix   = (size(xcen))[1]
   if NOT keyword_set(nycoeff) then nycoeff = 6
   if NOT keyword_set(ntcoeff) then ntcoeff = 4

   invvar = 1.0/xerr^2 * (xerr LT 90.0)
   y = ((2*dindgen(npix) - npix) / npix) # replicate(1,ntrace)

   if NOT keyword_set(t) then begin
     nrm0 = (ntrace-1)
     t = (2*dindgen(ntrace)/nrm0 - 1 ) ## replicate(1,npix)
   endif

   work2d = dblarr(npix*ntrace,nycoeff*ntcoeff)
;   worky = fchebyshev(y[*], nycoeff)
;   workt = fchebyshev(t[*], ntcoeff)
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

