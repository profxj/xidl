function bcextract, image, invvar, waveimg=waveimg, wset=wset, $
       tset=tset

   ncol = (size(image))[1]
   nrow = (size(image))[2]
   temparr = fltarr(ncol)

   temps = {     flux : temparr, $
             fluxivar : temparr, $
                 box  : temparr, $
             boxerr   : temparr, $
                 wave : temparr, $
              trace   : temparr, $
                width : temparr, $
           background : temparr }


    work = transpose(image)
    wivar = transpose(invvar)
    ximage = findgen(ncol) # replicate(1,nrow)
    yimage = replicate(1,ncol) # findgen(nrow)

;
;    try to set up wavelength scale.
;    
    if NOT keyword_set(waveimg) then begin
      if keyword_set(wset) then traceset2xy, wset, pixnorm, waveimg $
      else waveimg = ximage
    endif

    wave = transpose(waveimg)

;
;  Let's find the object
;
   sectorwidth = 100
   ns = ncol/sectorwidth

   medimage = djs_median(image,2) # replicate(1,nrow)
   bareimage = image - medimage
   sector = djs_median(reform(bareimage[0:ns * sectorwidth - 1,*], $
         sectorwidth,ns,nrow),1)

   print, max(sector,pl)
   xstart =  pl / ns
   ystart =  ((pl mod ns) + 0.5) * sectorwidth

   xfirst = trace_crude(work, xstart=xstart, ystart=start, yset=yfirst)
   xy2traceset, yfirst, xfirst, firstset, ncoeff=4, yfit=xfit1
   flux1 = extract_boxcar(transpose(bareimage),xfit1)

;
;  estimate the width of the strongest profile
;

   diff = yimage - xfit1 # replicate(1,nrow)
   close = where(abs(diff) LT 5.0)
   if close[0] EQ -1 then return, 0

   x = diff[close] 
   medimage = flux1 # replicate(1,nrow)
   y = (bareimage/medimage)[close]
   hmm = gaussfit(x,y,aa,nterms=3,estimates=[2.0,0.0,0.5])
   sigma = aa[2]

   
   fullset = { func : firstset.func, $
               xmin : firstset.xmin, $
               xmax : firstset.xmax, $
              coeff : firstset.coeff[*,0] # replicate(1,nrow) }        
   fullset.coeff[0,*] = findgen(nrow)

   pixfull=0
   traceset2xy, fullset, pixfull, fullcen

;
;  Now let's try to find any other apertures, to extract or mask out.
;
   boximage = extract_boxcar(transpose(bareimage), fullcen, radius=1.0)

;
;  Scrunch again
;

  boxsector = total(djs_median(reform(boximage[0:ns * $
         sectorwidth - 1,*], sectorwidth,ns,nrow),1),1)
  djs_iterstat, boxsector, median=box_m, sigma=box_s
  goodrows = (abs(boxsector - box_m) LT 3*box_s)  # replicate(1,ncol)
 
; third order fit:
    
  xy2traceset, findgen(nrow) # replicate(1,ncol) , $
        transpose(boximage), smooth_set, $
        mask=goodrows, yfit=smooth_fit, lower=2, upper=2

;  kernel = gauss_kernel(sigma,hpix=10*sigma) - $
;           gauss_kernel(3.0*sigma, hpix=10.0*sigma)

   radius = 3.0   ; 6 pixel boxcar ??
   temps.flux = extract_boxcar(transpose(bareimage)-smooth_fit,xfit1,$
        radius=radius)
   temps.wave = extract_boxcar(transpose(waveimg),xfit1,radius=radius)/6.
   return, temps
end
