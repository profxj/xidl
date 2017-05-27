pro guidefind, imagename, gg, fwhm=fwhm, $
      findsig=findsig, maxnstar=maxnstar, useatv=useatv, flux=flux, $
      outfwhm=outfwhm, mirror=mirror
;+
; NAME:
;   GUIDEFIND
;
; PURPOSE:
;   Find objects on guider image
;
; CATEGORY:
;   Spectroscopy
;
; CALLING SEQUENCE:
;   guidefind,image,gg, [findsig=findsig, 
;                            maxnstar=maxnstar,/USEATV,outfwhm=outfwhm,
;                             mirror=mirror]
;
; 
; INPUTS:
;   image -- either a 2d array containing a guider image OR the name
;            of a fits file containing the image
;
; KEYWORD PARAMETERS:
;   fwhm = fwhm -- the seeing fwhm in pix (Default:3); used to
;                  determine the smoothing kernel for object finding
;   findsig=findsig -- the number of sigma significance required for
;                      an object (Default: 5)
;   maxnstar=maxnstar -- if set, the max. # of objects desired (only
;                        the brightest will be returned)
;   /USEATV -- set to generate an ATV plot of the image with the
;              objects found marked (with different colors for the two
;              centroids - blue = CNTRD, red=ATV algorithm.  In
;                          general, the two are so close that only the
;                          red symbol is visible.
;   outfwhm=outfwhm -- returns estimated FWHM of the objects found, in
;                      pixels
;   mirror -- if set, keep only objects in mirror region
;
; OUTPUTS:
;   gg- structure containing:
;    xout,yout -- arrays containing x & y centroids (from FIND/CNTRD)
;                 for each object found
;    flux -- integrated flux returned as a keyword
;
; OPTIONAL OUTPUTS: (in gg also)
;    cxout,cyout -- arrays containing x & y centroids (from
;                   traditional centroid as implemented in ATV) for
;                   each object found
;
; MODIFICATION HISTORY:
;    11jun02 JAN
;    7jul02 MD - structure for output
;-



if n_elements(useatv) eq 0 then useatv = 0

if n_elements(findsig) eq 0 then findsig = 5 ; DAO find sigma cut
; following code from ptastrom.pro


s = size(imagename)

if s(0) ne 2 then begin
   image = readfits(imagename)
   s = size(image)
endif else image = imagename

nx = s[1]
ny = s[2]

maxsep = 10.*nx/1024. ; pixels



nofwhm = n_elements(fwhm) eq 0
if n_elements(fwhm) eq 0 then fwhm = ceil(7.*nx/1024.)
if nofwhm then print, 'FWHM used:', fwhm
on_error, 2

xarr = lindgen(nx, ny) MOD nx
yarr = lindgen(nx, ny) / nx

topbadslope = (489-500.)/498.*nx/512.
topbadint = 500.*nx/512.

midbadslope = (-10.)/495.*nx/512.
midbadint = 304.*nx/512.

badslopert = 32./11.*nx/512.
badintrt = (326-badslopert*511.)*nx/512

transy = (midbadslope*xarr+midbadint) > (badslopert*XARR+badintrt)

; define the two regions
bottom = where(yarr lt transy-10)
top = where(yarr gt transy+10)
bad =  where(yarr ge transy-10 AND yarr le  transy+10)
bad = [bad, where(yarr gt topbadslope*xarr+topbadint-5.)]

; create new image, with 0 median level and roughly constant gain
imnew = float(image)

; determine 'gain' levels on and off mask
pixforstatt = randomu(seed, 4001)*n_elements(top-1)
pixforstatb = randomu(seed, 4001)*n_elements(bottom-1)

medbottom = median(imnew[bottom[pixforstatb]])
medtop = median(imnew[top[pixforstatt]])

imnew(bottom) = (imnew(bottom)-medbottom)*medtop/medbottom
imnew(top) = imnew(top)-medtop
imnew = imnew-median(imnew, 7)
imnew[bad] = -0.001

; calculate statistics separately on top & bottom of camera
pixforstatt = randomu(seed, 40001)*n_elements(top-1)

djs_iterstat, imnew[top[pixforstatt]], sigrej=5, sigma=imsigt, mean=immeant, maxiter=3

pixforstatb = randomu(seed, 40001)*n_elements(bottom-1)

djs_iterstat, imnew[bottom[pixforstatb]], sigrej=5, sigma=imsigb, mean=immeanb, maxiter=3


; these parameters control DAOPHOT
; to select stars only:
;   sharplim = [.2, 1.0]
;   roundlim = [-1.0, 1.0]*1.5
;   hmin = imsig*findsig
;to select stars and galaxies:

   sharplim = [0.1, 1.1] ; allow galaxies, not cosmics
;   sharplim = [0., 1.4] ; tweak
   roundlim = [-1.0, 1.0]*2.0 ; allow elliptical objects

   imsig = imsigt > imsigb
   

   maskedpix = where(imnew lt -2.5*imsig, maskct)
   if maskct ne 0 then imnew[maskedpix] = -0.001
   crim = median(imnew, 3)
   whcr = where(imnew-median(imnew, 3) gt 5*imsig, crct)
   if crct gt 0 then imnew[whcr] = crim[whcr]
  
whtgood = where(imnew[top] ne -0.001)
whbgood = where(imnew[bottom] ne -0.001)

pixforstatt = whtgood[randomu(seed, 40001)*n_elements(whtgood-1)]


pixforstatb = whbgood[randomu(seed, 40001)*n_elements(whbgood-1)]


djs_iterstat, imnew[top[pixforstatt]], sigrej=5, sigma=imsigt, mean=immeant, maxiter=3


djs_iterstat, imnew[bottom[pixforstatb]], sigrej=5, sigma=imsigb, mean=immeanb, maxiter=3

   signoise = imnew
   signoise[top] = imnew[top]/imsigt
   signoise[bottom] = imnew[bottom]/imsigb

   hmin = findsig

   nobj = 1
   imx =[1, 2]
   
   fixfwhm:
   if (nobj eq 0 or max(imx) eq min(imx)) and fwhm lt 15./1024*nx  then begin
      fwhm = fwhm+1     
      message, 'increased fwhm to '+string(fwhm), /reset, /informational
   endif

   find, signoise, imx, imy, flux, sharp, round, hmin, fwhm, roundlim, $
     sharplim, /silent

   nobj = n_elements(imx)
   
   if nobj eq 0 or max(imx) eq min(imx) and fwhm lt 15./1024*nx then begin
      message, 'NO STARS YET', /informational
      goto, fixfwhm
   endif

   imcenx = imx
   imceny = imy

   bigbox = fix(1.5*fwhm+0.5)
   cenbox = fix(1.8*fwhm+0.5)

   fwhmarr = fltarr(nobj)

   for i=0, nobj-1 do begin
      split_radplotf, imnew, imnew eq -0.001, imx[i], imy[i], tmpfwhm
      fwhmarr[i] = tmpfwhm
      imcenterjan, imnew, imx[i], imy[i], xcen, ycen, big=bigbox, cbox=cenbox
      wayoff = (abs(xcen-imx[i]) > abs(ycen-imy[i])) gt (bigbox*0.7)
      if NOT wayoff then imcenx[i] = xcen
      if NOT wayoff then imceny[i] = ycen
;      if wayoff then print, 'centroid failure - set to match CNTRD'
   endfor


   if keyword_set(useatv) then if useatv gt 0 then begin
      atv, imnew
      atvplot, imx, imy, psym=1, color=3
      atvplot, imcenx, imceny, psym=1, color=1
   endif


; the below error-checking code comes from ptastrom



; when there are no good stars returned, an array of zero flux stars
; may be returned (I think this is a DAO find bug).  -DPF
   count= 0
   if n_elements(flux) GE 1 then goodstar = where(flux GT 0, count)

;   IF count LT 5 THEN BEGIN 
;      errstr = string('---> DAO Find:', count, ' stars found -- giving up!', $
;                      format='(A,I4,A)')
;      pt_errlog, infostr, errstr
;      return
;   ENDIF 

; check if too many stars
   IF keyword_set(maxnstar) THEN BEGIN 
      nstar = n_elements(imx)
      IF (nstar GT maxnstar) THEN BEGIN 
         ncull = (nstar/4) > maxnstar    ; never use fewer than 1/4 of stars
         print, 'FIND: found', nstar,' stars - culling to', ncull
;         ind = long(findgen(ncull)*float(nstar)/float(ncull))
         sind = sort(flux)
         ind = sind[nstar-ncull-1:nstar-1] ; keep brightest ncull
         imx = imx[ind]
         imy = imy[ind]
         flux = flux[ind]
         sharp = 0  ; throw away
         round = 0  ; throw away
      END
   END 

   ggt = {xout: 0., yout: 0., xsout: 0.,  ysout: 0., $
           flux: 0.} ;assign structure for returning results
   gg = replicate(ggt, n_elements(imx))
   gg.xout = imx
   gg.yout = imy
   gg.xsout = imcenx
   gg.ysout = imceny
   gg.flux =  flux


   outfwhm =  fwhmarr 
return
end




