

function mike_fweight, img, xcen, ycen, invvar=invvar, $
                       niter=niter, radius=radius, xerr=xerr, sig=sig

     if NOT keyword_set(sig) then sig=2.0

     sz_xcen = size(xcen,/dimen) 

     if NOT keyword_set(ycen) then begin
        ycen = n_elements(sz_xcen) EQ 1 ? findgen(sz_xcen[0]) : $
                findgen(sz_xcen[0]) # replicate(1,sz_xcen[1])
     endif

     if NOT keyword_set(niter) then niter=1L

     xtmp = xcen
     for j=1,niter do xtmp = trace_fweight(img, xtmp, ycen, $
                             radius=radius, xerr=xerr_fweight)

     ;
     ;  Iterations are done, estimate errors based on flux and variance
     ;
     fluxmid    = extract_boxcar(img, xtmp, radius=radius)
     fluxleft   = extract_boxcar(img, xtmp-radius, radius=radius)
     fluxright  = extract_boxcar(img, xtmp+radius, radius=radius)
     badrow = xcen*0.
     if keyword_set(invvar) then begin

       var = 1/(invvar + (invvar EQ 0)) * (invvar GT 0)
       badmask = 1.0*(invvar LE 0)
       badrow =  extract_boxcar(badmask, xtmp, radius=radius)
       varleft  = extract_boxcar(var, xtmp-radius, radius=radius)
       varright = extract_boxcar(var, xtmp-radius, radius=radius)
       var_numer =  4.0*(varleft>varright)  
       var_numer = var_numer + (var_numer EQ 0) 
     endif else var_numer =  4.*(fluxleft > (fluxright > 0)) + 10.
 

     err_denom =  (2.0*fluxmid - fluxleft - fluxright) 
     inv_err = (err_denom > 0) /sqrt(var_numer)
     xerr = 1.0 / (inv_err > 1./radius) 

     for iord=0,n_elements(sz_xcen) EQ 1 ? 0 : sz_xcen[1]-1 do begin
        djs_iterstat, err_denom[*,iord], mask=om, sigrej=5
        badrow[*,iord] = badrow[*,iord] + (om EQ 0)
     endfor

     baderr= where(inv_err LE 1./radius OR (badrow GT 0) OR $
                   (err_denom LT sig* sqrt(var_numer)))
     if baderr[0] NE -1 then xerr[baderr] = 999
    
     return, xtmp 
end
     
