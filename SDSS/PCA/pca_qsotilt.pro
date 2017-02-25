;+ 
; NAME:
; pca_qsotilt   
;   Version 1.0
;
; PURPOSE:
;  Attempt to fit in logarithmic flux a low order polynomial to
;   match the CIV region (range)
;   coeffs contains the [npoly,nspec] coefficients based on rlam
;
;   This fit is done in logarithmic flux vs. logarithmic wavelength
;
; CALLING SEQUENCE:
; pca_qsotilt, rlam, flux, invvar, xmed, ymed, tilt, tiltivar, $
;                coeff, npoly=, range=, answer=, FLG_TILT=
;
; INPUTS:
;  rlam -- Observed Wavelength array (logarithmic)
;  flux -- Flux of quasar
;  invvar -- Inverse variance array
;  xmed -- Rest frame wavelength array [logarithmic; rlam / (1+zem) ]
;  ymed -- 0th channel of the PCA fit  (pca[*,0])
;
;
; RETURNS:
;  tilt --  Normalized fit to QSO (tilt taken out)
;  tiltivar -- Inverse variance of the fit
;  coeff -- Coefficients of the fit to the tilt
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  npoly=  -- Number of coefficients in POLY fit to tilt
;  range=  -- Rest-wavelength range to fit [default: 1270--2000]
;  pixlim= -- Minimum pixels necessary for good fit [default: 10]
;
; OPTIONAL OUTPUTS:
;  ANSWER= -- Tilted fit
;  FLG_TILT -- Flag describing the success of the fit (0=bad, 1=good)
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Sep-2003 Written by JXP (based on qso_tilt by SB)
;-
;------------------------------------------------------------------------------
pro pca_qsotilt, rlam, flux, invvar, xmed, ymed, tilt, tiltivar, $
                 coeff, npoly=npoly, range=range, answer=answer, $
                 FLG_TILT=flg_tilt, PIXLIM=pixlim

  if  N_params() LT 6  then begin 
    print,'Syntax - ' + $
             'pca_qsotilt, rlam, flux, invvar, xmed, ymed, tilt, tiltivar, coeff, NPOLY='
    print, '    RANGE=, ANSWER=, FLG_TILT=  [v1.1]'
    return
  endif 


;   nspec = (size(flux))[2] 
   npix = n_elements(rlam)
   tilt = fltarr(npix)
   tiltivar = fltarr(npix)

   if NOT keyword_set(npoly) then npoly=2
   if NOT keyword_set(range) then range = alog10([1270., 2000.])
   if not keyword_set(pixlim) then pixlim = 10
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
       if ct2 GT pixlim then begin
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
        
        




   
   
      
