;+ 
; NAME:
; x_findgauss   
;   Version 1.0
;
; PURPOSE:
;    Fits a continuum to spectroscopic data interactively
;
; CALLING SEQUENCE:
;   
;   dla_sdssrich, fil
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   parse_sdss, fil
;
;
; PROCEDURES/FUNCTIONS CALLED:
; REVISION HISTORY:
;   24-Sep-2003 Written by SB
;-
;------------------------------------------------------------------------------
;  Attempt to fit in logarithmic flux a low order polynomial to
;   match the CIV region (range)
;   coeffs contains the [npoly,nspec] coefficients based on rlam
;
;   This fit is done in logarithmic flux vs. logarithmic wavelength


pro x_findgauss, y, invvar, xarr=xarr, sigma=sigma, $
   xpeak=xpeak, ypeak=ypeak, xsn=xsn, sn=sn, $
   width=width, nfind=nfind, minsep=minsep


   if NOT keyword_set(sigma) then sigma = 1.0
   if NOT keyword_set(width) then width = long(2.5*sigma)
   if NOT keyword_set(minsep) then minsep = long(4.1*sigma)
   if NOT keyword_set(nfind) then nfind = 50L
  
   if NOT keyword_set(invvar)  then invvar = y*0 + 1.
   if NOT keyword_set(xarr)  then  xarr = lindgen(n_elements(y))

   gkernel = gauss_kernel(sigma)
   gaussivar = convol(invvar, gkernel^2, /edge_truncate)
   gnumerator = convol(y*invvar, gkernel, /edge_truncate)

   gfilter = gnumerator / (gaussivar + (gaussivar LE 0)) * (gaussivar GT 0)
   gsn = gfilter * sqrt(gaussivar) 

   xsn = find_npeaks(gsn, xarr, nfind=nfind, $
                 width=width, minsep=minsep, ypeak=sn) 

   ;; Avoid zero spots
;   print, n_elements(xsn)
   gd = where(invvar[round(xsn)] GT 0., ngd)
   if ngd EQ 0 then xpeak = 0L else xsn = xsn[gd]

   ;; Tune
   xpeak = xsn *0.0
   ypeak = xsn*0.0
   for j=0,n_elements(xpeak)-1L do begin 
      xpeak[j] = zfitmax(gfilter, xguess=xsn[j], ypeak=tempy, width=width) 
      ypeak[j] = tempy 
   endfor
         

   return
end

   
