;+ 
; NAME:
; x_findgauss   
;   Version 1.0
;
; PURPOSE:
;    Finds Gaussian absorption features in a spectrum (should be
;    normalized).  The routine uses find_npeaks to find the features
;    and then tunes up on then with zfitmax.
;
; CALLING SEQUENCE:
;  x_findgauss, y, [invvar], xarr=, sigma=, xpeak=, ypeak=, xsn=, sn=, 
;  width=, nfind=, minsep=
; 
; INPUTS:
;  y  -- Flux
;  [invvar] -- Inverse variance
;
; RETURNS:
;
; OUTPUTS:
;  xpeak -- x values of the peaks
;  ypeak -- y height of the peak
;
; OPTIONAL KEYWORDS:
;  SIGMA= -- Width of Gaussian (pixels)  [default: 1.]
;  WIDTH= -- FWHM of the feature [default: 2.5*sigam]
;  MINSEP= -- Minimum separation between lines [default: 4.1*sigma]
;  NFIND=  -- Number of peaks to find with find_nepaks [default: 50]
;  XARR=  -- x array (e.g. wavelength)
;
; OPTIONAL OUTPUTS:
;  xsn -- x values of peaks from find_npeaks
;   sn -- y values of the peaks from find_npeaks
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
; find_npeaks
; gauss_kernel
; zfitmax
;
; REVISION HISTORY:
;   24-Sep-2003 Written by SB
;   2003  Grabbed by JXP and slightly modified
;   15 Jul 2015 Added /raw option, KLC
;-
;------------------------------------------------------------------------------

pro x_findgauss, y, invvar, xarr=xarr, sigma=sigma, sigxpeak=sigxpeak, $
                 xpeak=xpeak, ypeak=ypeak, xsn=xsn, sn=sn, gsn=gsn, $
                 width=width, nfind=nfind, minsep=minsep, raw=raw

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
          'x_findgauss, y, [invvar], xarr=, sigma=, xpeak=, ypeak=, xsn=, sn=, width=, nfind=, minsep= [v1.1]'
    return
  endif 

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

   xsn = x_findnpeaks(gsn, xarr, nfind=nfind, $
                      width=width, minsep=minsep, ypeak=sn) 


   ;; Avoid zero spots
   gd = where(invvar[round(xsn)] GT 0., ngd)
   if ngd EQ 0 then xpeak = 0L else xsn = xsn[gd]

   if keyword_set(raw) then begin
      ;; For low-res, low-SNR (or just desiring less computationally
      ;; intensive) return the raw results
      if ngd eq 0 then begin
         sigxpeak = 0.
         ypeak = 0.
      endif else begin
         xpeak = xsn
         sigxpeak = replicate(width,ngd)
         ypeak = gfilter[xpeak]
      endelse
      
   endif else begin
      ;; Tune
      xpeak = xsn *0.0
      sigxpeak = xpeak
      ypeak = xsn*0.0
      for j=0,n_elements(xpeak)-1L do begin 
         xpeak[j] = zfitmax(gfilter, xguess=xsn[j], ypeak=tempy, width=width,xerr=tempxerr)
         sigxpeak[j] = tempxerr
         ypeak[j] = tempy 
      endfor
   endelse

   return
end

   
