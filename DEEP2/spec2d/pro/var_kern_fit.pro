;+
; NAME:
; var_kern_fit
;
; PURPOSE:
; Finds kernel to convolve a spectrum f1 to f2 in local windows, then
; fits the kernel to a polynomial across the domain and returns the
; convolved spectrum f1prime.
;
; CALLING SEQUENCE:
; 
; f1prime = var_kern_fit(x, f1, f2, measerr, window, nhalf, npoly)
;
; INPUTS:
;
; x       - independent variable
; f1      - the unconvolved template spectrum to be convolved
; f2      - the spectrum to which we want to convolve f1
; measerr - measurement errors in f2
; window  - size of window in which to "measure" the local kernel
; nhalf   - half width of kernel (the kernel is an array of length
;           [2*nhalf]+1). 
; npoly   - The (order+1) of the polynomial used to fit the local
;           "measurements" of the kernel as a function of x.
;
; OPTIONAL INPUTS:
;
; threshold - threshold above which a signal is considered significant
;             in f1 (used for deweighting continuum regions).  Default
;             threshold is 1/50 of the max in the local window, or a
;             6-sigma deviation from the continuum, whichever is
;             greater.  Set to -1 to use the whole spectrum equally
;             (or set deweight to 1).
; contin1   - Use this to feed the routine an explicit model for the
;             continuum of f1 (default: median smooth over 2*window
;             box).
; contin2   - same, as above, except for f2.
; contin2prime - continuum model for convolved function.  Default=contin2
; conterr2  - error in contin2; default is Poisson errors.
; deweight  - Factor by which to increase the variance in
;             continuum regions 
; winstep   - The amount by which the fit window is stepped while
;             finding kernels as a function of wavelength.  Default is
;             window/2. 
;	
; KEYWORDS:
;
; contsub    - if set, continuum models are subtracted from f1 and f2
;             before fitting; non-subtracted funtions are returned.
; medkern   - if set, use median value of kernel, instead of
;             polynomial fit.
; spline    - if set, use spline interpolation to smooth kernel,
;             rather than polynomial fitting
;
; OUTPUTS:
;
; f1prime   - f1 convolved with variable kernel.
;
; OPTIONAL OUTPUTS:
;
; fullsignal- Array indicating at what x values significant signal was
;             found (zero for continuum, nonzero for signal).
; kernels   - Measurements of the kernel elements for each of the nfit
;             windows used in the fit.  Dimensions: nfit x 2*nhalf+1
; xfit      - Central x value of each window used in the fit
; fitkern   - The polynomial fit for each kernel element, evaluated at
;             xfit.  Dimensions:  nfit x 2*nhalf+1
; fft       - if set, uses FFT technique to find kernel
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
; Written by BFG summer 2002
; Made presentable 14Oct02
; Revised by JAN Nov 02
; Revised by BFG 11Nov2002
;
;----------------------------------------------------------------------

function var_kern_fit, x, f1, f2, measerr, window, nhalf, npoly, $
            winstep=winstep, threshold=threshold, contin1=contin1, $
            contin2=contin2, cont2prime=cont2prime, fft=fft, $
            conterr2=conterr2, deweight=deweight,  $
            fullsignal=fullsignal,  contsub=contsub, $
            medkern=medkern,  kernels=kernels, $
            xfit=xfit,  fitkern=fitkern,  spline=spline


;---find domain
  domain = minmax(x)

;---make continuum models if necessary

  if not keyword_set(contin1) then $
    contin1 = djs_median(f1, width=2*window, boundary='reflect')
  if not keyword_set(contin2) then $
    contin2 = djs_median(f2, width=2*window, boundary='reflect')
  if not keyword_set(cont2prime) then cont2prime = contin2

  if not keyword_set(conterr2) then $
    conterr2=sqrt(contin2)

;---definitions
  if not keyword_set(winstep) then $
    winstep = window/2.
  nregions =  (domain[1]-domain[0])/winstep ; number of fit regions
                                ;NB: these are overlapping tiles of length window.
  fullsignal =  intarr(n_elements(x)) ; where there's significant signal
  kernels =  fltarr(nregions, 1+2*nhalf) ; kernel "measurements"
  sigout =  kernels             ; error in kernel "measurements"
  dofit = intarr(nregions)      ; where in x to do the fit
  xfit =  fltarr(nregions)      ; x-values for kernel "measurements"
  if NOT keyword_set(spline) then $
    kfit =  fltarr(1+2*nhalf,  npoly)  $; array to store polynomial coeff.
  else $                        ;alternatively, store spline coefficients
    kfit=   fltarr(1+2*nhalf, nregions)
  fitkern =  fltarr(nregions, 1+2*nhalf) ; fit kernel in each region

;---find regions with signal; broaden slightly to get wings
     if not keyword_set(threshold) then begin $
       threshold =  max(f1)/50.+contin1 > contin1+6.*sqrt(contin1)
     endif else if threshold[0] eq -1 then begin
        signal =  fltarr(n_elements(f1))+1
     endif

     if threshold[0] ne -1 then begin
        fullsignal =  (f1 gt threshold)
        
        if total(fullsignal) lt 20*nhalf then begin
            threshold=threshold-15
            message,'threshold was lowered!!',/INFO
            fullsignal =  (f1 gt threshold)
        endif
        fullsignal=fullsignal OR dilate(fullsignal,intarr(3*nhalf)+1)

     endif

  l1 = n_elements(f1)
;---first set up weighting array to apodize spectra (10% at each end)
    h = hanning(4*nhalf+1)
     weight = fltarr(l1)+1
     weight[0:2*nhalf] = h[0:2*nhalf]
     weight[l1-2*nhalf-2:l1-2] = h[2*nhalf:4*nhalf]
     weight[l1-1]=0.


;---apodize spectra
  fit1 = (f1-contin1)*weight
;  fit2 = (f2-contin2)*weight
;---set up array of shifts
  shifts = fltarr(2*nhalf+1, l1)
  for i=-nhalf, nhalf do begin
     shifts(i+nhalf, *) = shift(fit1, -i)  
    ;to keep the kernel defined as for the fft, we need this choice of sign
  endfor   


;---march across domain, overlapping window.
  
  for i=0, nregions-1 do begin

;---choose region to convolve
     xmin =  domain[0]+i*winstep
     xmax =  xmin+window
     if xmax ge domain[1] then begin 
        xmin =  domain[1]-window
        xmax = domain[1] 
     endif
     xfit[i] =  (xmax+xmin)/2.  ; x value for polynomial fit
     ireg =  where((x ge xmin) and (x le xmax))
     xreg = x[ireg]
     f2reg =  f2[ireg]
     f1reg =  f1[ireg]
     shiftreg = shifts[*, ireg]
     err =  measerr[ireg]
     signal = fullsignal[ireg]

;---subtract off continuum from both functions (must do this here so
;   that threshold is set correctly above).
     if keyword_set(contsub) then begin
        f1reg = f1reg-contin1[ireg]
        f2reg = f2reg-contin2[ireg]
     endif

;---deweight continuum regions by increasing errors

     if not keyword_set(deweight) then deweight = 1.

     icont =  where(signal eq 0,contct)
     if contct gt 0 then err[icont] = sqrt(deweight)*err[icont]
     measerr[ireg] =  err


;---if there are enough pixels with signal to get a reasonable fit,
;   find kernel 

     if total(signal) ge (6*nhalf) then begin
        dofit[i] =  1
;        print,  xfit[i], dofit[i]
;---find kernel to go from sky to object spectrum in each region          
        kernels[i, *]=  find_kernel_reg(shiftreg,  f2reg,  nhalf, $
                             err, fit,  sigouti,  const,  fft=fft) 
        sigout[i, *] =  sigouti ;error in the fit for this window
     endif

;  print, const
  endfor

  wfit =  where(dofit,fitct)  
  if fitct eq 0 then begin
      return,-1
  endif

  if keyword_set(medkern) AND fitct gt 1 then begin
;---set kernel elements to the median of the various measurements
     for i=0, 2*nhalf do begin
        kernelsi =   reform(kernels[wfit, i])

        djs_iterstat, kernelsi,  median=medkerni
        kernels[wfit, i] = medkerni
     endfor
     
  endif

  mskerns = kernels ;to store median smooth below
  sumkernel = total(kernels, 2) ;total of kernel for each lambda
  sdkerns = fltarr(2*nhalf+1);to store stddev of kernels
;---fit each pixel of kernel to a polynomial in the lambda direction
  for i=0, 2*nhalf do begin  
;---first, reject outliers with a median smooth
     mskerns[*, i] = djs_median(kernels[*, i], width=3,boundary='reflect')
     djs_iterstat, kernels[*, i], sigma=sdk
     sdkerns[i] = sdk
     for j=0, n_elements(wfit)-1 do begin
        if kernels[j, i] gt mskerns[j, i]+6*sdkerns[i] then $
           sigout[j, i]=100*sigout[j, i]
     endfor

     sigout[0,*]=10.*sigout[0,*]
     sigout[n_elements(wfit)-1,*]=10.*sigout[n_elements(wfit)-1,*]

     kernelsi =   reform(kernels[wfit, i])

;---do polynomial fit.

;---first, remap domain to [-1,1]
     xx = (xfit[wfit]-domain[0])/((domain[1]-domain[0])/2.) - 1
     spl_x = xfit + findgen(n_elements(xfit))*.01 ;something to keep spline initiation happy--a big kludge for the moment, but x values must be ascending in the spline process
;     print,  xx

     if NOT keyword_set(spline) then begin ;do polynomial fit
       kfit[i, *] =  svdfit(xx,  kernelsi, npoly, $
                            measure_errors=sigout[wfit, i], $
                            yfit=fitkerni, $
                            /double, /legendre) ;double is necessary.
       fitkern[wfit, i] =  fitkerni
     endif else $ ;spline fit
       kfit[i, *] = spl_init(spl_x, kernelsi) ;initialize spline
  endfor
  
;---subtract continua from entire spectra (have to do this after so
;   that threshold is set correctly above)
  if keyword_set(contsub) then begin
     f2 = f2-contin2
     f1 = f1-contin1
  endif

;---convolve to get fit.
  if NOT keyword_set(spline) then begin ;polynomial method
    if NOT keyword_set(fft) then $  ;if FFT method, use prewhitening
       f1prime =  convol_var(x,f1, nhalf,  npoly, kfit) $
    else $
       f1prime =  convol_var(x,f1/measerr, nhalf,  npoly, kfit)*measerr
  endif else begin  ;spline method
    if NOT keyword_set(fft) then $  ;if FFT method, use prewhitening
       f1prime =  convol_spl(x,f1, spl_x, kernels,  kfit) $
    else $
       f1prime =  convol_spl(x,f1/measerr, spl_x, kernels , kfit)*measerr
  endelse
;---add medians back
  if keyword_set(contsub) then begin
     f2 =  f2+contin2
     f1 =  f1+contin1
     f1prime = f1prime+cont2prime  
  endif


;---only return values used in the fit
  xfit = xfit[wfit]
  kernels = kernels[wfit, *]
  fitkern = fitkern[wfit, *]

  return, f1prime

end









  
