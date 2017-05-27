;+
; NAME: find_kernel_reg
;
; PURPOSE: Performs a linear regression fit to find the discrete
; kernel needed to convolve one function (f1) into another (f2).
;
; CALLING SEQUENCE:
; cc = find_kernel_reg(shiftarr, f2, nhalf, measerr, fit, sigout, $
;                const, fft=fft)
; 
; INPUTS:
; 
; shiftarr- A (2*nhalf+1 x n_elements(f2)) array giving the function to
;           be convolved, shifted by nhalf in either direction 
; f2      - The desired result of the convolution (f2 = cc * f1, where
;           * is the convolution symbol)
; nhalf   - the half-width of the kernel to be returned (the kernel is
;           an array of length 2*nhalf + 1)
; measerr - the measurement errors in the function f2
; fit     - the result of convolving f1 with the returned kernel cc
;           (should approximate f2 well)
; sigout  - estimate of the error in each element of the kernel.
; const   - constant term returned from the regression
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;   fft   - set on if FFT technique desired for solution
;
; OUTPUTS:
;
; cc      - The kernel to convolve between f1 and f2 (an array of
;           length [2*nhalf]+1)
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
; 
; Written by JAN back in the mists of time.
; Made presentable by BFG on 14Oct02
; FFT method by MD 11Nov02
;
;----------------------------------------------------------------------

function find_kernel_reg, shiftarr, f2, nhalf,measerr, fit, sigout, $
         const,  fft=fft


; do a regression to fit this array to the function f2 with
; coefficients given by cc (i.e., find the kernel) 
; DOUBLE IS NECESSARY!
  if NOT keyword_set(fft) then begin
    cc = regress(shiftarr, f2,  meas=measerr, sigma=sigout, yfit=fit,$ 
               const=const, /double)

    cc = (reform(cc, n_elements(cc)))

  endif else begin ;do FFT technique instead

    ssize = n_elements(f2)
    nshifts = (size(shiftarr, /dimen))[0]
    nsky = reform(shiftarr[nshifts/2, *]) 
;    tnoise = total(measerr^2) ;accumulate total noise power
;    pkerr = sqrt(tnoise/(ssize)) ;noise per wavenumber
    pkerr = 1. ;noise in normalized case should be 1/mode
; test reducing noise
; make noise properties more uniform by dividing by measerr
    nskyn = nsky/measerr
    f2n = f2/measerr

;extract central example, which is unshifted

    h = hanning(ssize/5)
    nhan = n_elements(h)/2 ;half window
    f2e = [f2n, fltarr(ssize)] ; pad array with zeros
    f2e[0:nhan-1] = f2e[0:nhan-1]*h[0:nhan-1] ;hanning smooth
    f2e[ssize-nhan:ssize-1] = f2e[ssize-nhan:ssize-1]*h[nhan:2*nhan-1]
    nskye = [nskyn, fltarr(ssize)] ;pad skyslit as well
    nskye[0:nhan-1] = nskye[0:nhan-1]*h[0:nhan-1] ;hanning smooth
    nskye[ssize-nhan:ssize-1] = nskye[ssize-nhan:ssize-1]* $
               h[nhan:2*nhan-1]
;hanning smooth outer 10% of each spectrum, each end
    f2ek = fft(f2e)
    nskyek = fft(nskye) ;get FFT of each array

; use a Weiner filter type argument in doing the division
    nskyxy = fltarr(2, n_elements(nskyek)) ;make array of x,y
    f2xy = nskyxy
    nskyxy[0, *] = float(nskyek) ;populate with real, imaginary
    nskyxy[1, *] = imaginary(nskyek)
    f2xy[0, *] = float(f2ek)
    f2xy[1, *] = imaginary(f2ek) 
    nskypolar = cv_coord(from_rect=nskyxy,  /to_polar) ;get amp. and phase
    f2polar  =  cv_coord(from_rect=f2xy,    /to_polar)

; get mean value for use in filter
    ampsky = djs_median(reform(nskypolar[1, *]), width=ssize/20, $
        boundary='reflect')
;    filter = sqrt(ampsky^2/(pkerr^2 +ampsky^2)) ;Weiner filter    
    filter = (ampsky^2/(pkerr^2 +ampsky^2))^(.25)  ;Weiner filter

    filter = filter/mean(filter[0:ssize/20]) ;normalize to unity at small k
  
   ratiop = f2polar*0. 
;do Weiner filter on smoothed noise 
    ratiop[1, *]= f2polar[1, *]/ nskypolar[1, *]* filter ;Weiner suppression corrected
; in the mean
    ratiop[0, *] = f2polar[0, *] - nskypolar[0, *] ;angles subtract
    ratio = cv_coord(from_polar=ratiop, /to_rect)
    ratio = reform(complex(ratio[0, *], ratio[1, *])) ;convert back to complex
;    ratio = djs_median(ratio, width=5, boundary='reflect')

    kernel = shift(fft(ratio, /inverse), ssize)/n_elements(ratio)
;inverse transform and shift
    cc = float(kernel[-nhalf+ssize:nhalf+ssize]) ;extract kernel
    dd = regress(shiftarr, f2,  meas=measerr, sigma=sigout, yfit=fit,$ 
               const=const, /double)
; do the regression to get the measurment error
;    sigout = cc*0. +1. ;equal errors per point for this fit

;    save, file='nl_test.sav' ;save for testing
;stop

  endelse

return, cc
end







