;+ 
; NAME:
; pca_qsotilt   
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
;   09-Sep-2003 Written by JXP (based on qso_tilt by SB)
;-
;------------------------------------------------------------------------------
;  Attempt to fit in logarithmic flux a low order polynomial to
;   match the CIV region (range)
;   coeffs contains the [npoly,nspec] coefficients based on rlam
;
;   This fit is done in logarithmic flux vs. logarithmic wavelength

pro pca_qsotilt, rlam, flux, invvar, xmed, ymed, tilt, tiltivar, $
                 coeff, diff, npoly=npoly, range=range, answer=answer, $
                 FLG_TILT=flg_tilt

;   nspec = (size(flux))[2] 
   npix = n_elements(rlam)
   tilt = fltarr(npix)
   tiltivar = fltarr(npix)

   if NOT keyword_set(npoly) then npoly=2
   if NOT keyword_set(range) then range = alog10([1270., 2000.])
;   coeff = fltarr(npoly)
;   diff  = fltarr(nspec)

   bottom = max([min(range),min(xmed)])
   top    = min([max(range),max(xmed)])

   tt = bspline_iterfit(xmed, ymed, everyn=20, yfit=yfit, /silent)

;   print, 'Tilting with a base range:', 10^range
      
; FIT
   good = where(rlam GE bottom AND rlam LE top $
                AND invvar GT 0.0,ct)
   if ct GT 10 then begin
       model = bspline_valu(rlam[good], tt)
       goodpix = where(flux[good] GT 0 AND model GT 0,ct2)
       if ct2 GT 10 then begin
           use = good[goodpix]
           coeff = func_fit(rlam[use], $
                                 alog10(flux[use]/model[goodpix]), $
                                 npoly, func='poly')
           answer = 10^poly(rlam[*], coeff)
           tilt[*] = flux[*] / answer
           tiltivar[*] = invvar[*] * answer^2
           flg_tilt = 1
       endif else stop
   endif else begin
       answer = replicate(1., npix)
       tilt[*] = flux[*] / answer
       tiltivar[*] = invvar[*] * answer^2
       flg_tilt = 0
   endelse

 ; Burles counter of column number...
;   print, format='($, ".",i4.4,a5)', i, string([8b,8b,8b,8b,8b])
   
   return
end     
        
        




   
   
      
