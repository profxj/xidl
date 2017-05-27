;+
; NAME: 
; convol_var
;
; PURPOSE:
; Convolves a (discrete) function f(x) into a function fprime(y) with a
; (discrete)kernel K(y-x) that varies as a polynomial function of x.
;
; CALLING SEQUENCE:
; fprime = convol_var(x, f, nhalf, npoly, kfit, medkern=medkern, $
;                     kernels=kernels) 
;
; INPUTS:
; 
; x      - The independent variable for f
; f      - The function to be convolved.
; nhalf  - The half-width of the kernel (the kernel is an array of
;          length [2*nhalf]+1) 
; npoly  - the (order+1) of the polynomial describing the kernel
;          elements' variation.
; kfit   - 2d array containing the polynomial coefficients for each
;          kernel element (dimensions: 2*nhalf+1 x npoly)
;
; OPTIONAL INPUTS:
; kernels - Array of N "measurements" of the kernel elements, whose
;           median is to be computed. Ignored if the
;           medkern keyword is not set. Dimensions: N x 2*nhalf+1.
;
; KEYWORDS:
; medkern - if set, use median values of kernel elements, rather than
;            poly fit.  Must also set kernels if this is set.
; spline  - if set, use spline interpolation rather than polynomial smoothing
;
; OUTPUTS:
; fprime  - The result of convolving f with the variable kernel.
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
; Written by BFG, summer 2002. 
; Made presentable, 14Oct02.
; Revised by BFG 11Nov02
;
;----------------------------------------------------------------------

function convol_var,  x, f,  nhalf, npoly, kfit,xx=xx, $
                      medkern=medkern,  kernels=kernels,  spline=spline


  npts = n_elements(x)  

  kernel =  fltarr(npts, 2*nhalf+1);to store kernel evaluated at each point x
  fprime = fltarr(npts); to store result
  for i=0, 2*nhalf do begin
    if keyword_set(medkern) then begin ; use median of kernel "measurements"
      if not keyword_set(kernels) then $
        print, 'Error: Need to set kernels in convol_var when medkern is set.'
      kernel[*, i] =  djs_median(kernels[*, i])
    endif else begin
     ; for j=0, npoly-1 do begin
       ; evaluate polynomial fit at each point
        xx = (x-x[0])/((x[npts-1]-x[0])/2.)-1 ;done in calling argument
        if NOT keyword_set(spline) then $
           kernel[*, i] = polyleg(xx, kfit[i, *]) $ ;polynomial expansion
        else $ 
           kernel[*, i] = spl_interp(xx, f, kfit[i, *], xx) ;not going to work
      ;endfor
    endelse
    fprime = fprime+(kernel[*, i] * shift(f, nhalf-i)) 
                    ;add the next term in the convolution
  endfor

return,  fprime
end




