;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; measureSN.pro
; Author: Kathy Cooksey     Date: 22 Feb 2005
; Project: OVI survey and HST Archive Metal-Line System survey with 
;          Jason Prochaska
; Description: Function to return most common flux and error 
;              (and signal-to-noise ratio from it) 
; Input: 
;   spec -- name of spectrum FITS or flux array (inflg = 4)
; Optional:
;   err -- name of error FITS or error array (inflg = 4)
;   wave -- wavelenth array (inflg = 4)
;   inflg -- type of FITS file format
;         0  STIS data (*_f.fits and *_e.fits)
;         1  Flux and error in same file
;         2  Flux, error and wavelength in same file
;         3  FUSE format
;         4  directly input arrays (err and wave must be set)
;   /plot -- calls x_splot to show histogram of flux and error
;            or continuum
;   /cont -- measure S/N on rough continuum
;   /silent -- don't print messages
; Output: 
;   [mdflux,mderror,snr] -- array of most common flux and error
;                         and S/N thereof
; Optional:
;   region=array -- return the indices of regions used to define
;                   continuum (applies only if /cont set)
;   dev -- deviation of output SNR array
; History:
;  22 Feb 05 -- created by KLC (tho not completed)
;   9 Sep 05 -- turned into function and completed
;  13 Sep 05 -- limit to flux > 0. and corresponding error
;  14 Sep 05 -- prevent histograms from being too large;
;               check for NaN problems; add /silent
;  21 Oct 05 -- use /cont to measure S/N on continuum
;  24 Oct 05 -- adde region option to return continuum
;  24 Feb 05 -- actually use autocontfit for /cont option
;  20 Jun 07 -- added dev function
;  26 Nov 07 -- buffer histograms if not enough to gaussfit
;   6 Apr 08 -- Enable input of flux, error, wave arrays
;  28 Apr 08 -- Don't mess with input arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function measureSN,spec,err=err,inflg=inflg,plot=plot,wave=wave,$
                   cont=cont,dev=dev,region=region,silent=silent
if not keyword_set(inflg) then inflg = 0
case inflg of
   4: begin
      flux = spec
      if not keyword_set(err) then stop,'measureSN: must set err'
      error = err
      if not keyword_set(wave) then stop,'measureSN: must set wave'
      svspec = {wave:wave,flux:spec,error:err} ; save to prevent mistakes
      npix = n_elements(wave)
   end
   else: flux = x_readspec(spec,inflg=inflg,head=hd,sig=error,wav=wave,$
                           fil_sig=err,npix=npix)
endcase 

;;For those freak cases where some element is NaN, eliminate it in
;;both flux and error
bd = where(finite(flux,/nan) or finite(error,/nan) or $
           error eq 0.,nbd,complement=gd)
if nbd ne 0 then begin
    flux = flux[gd]
    error = error[gd]
    wave = wave[gd]
    if not keyword_set(silent) then begin
       if size(spec,/type) eq 7 then $
          print,'measureSN: elements of flux or error NaN, trim ',spec $
       else print,'measureSN: elements of flux or error NaN, trim'
    endif                       ; /silent
endif 

if not keyword_set(cont) then begin
;;Histogram data and determine peak
    div = 5.

;;To prevent trying to histogram large range of flux and/or error
;;limit the number of histogram bins by setting sztol
    sztol = 100*n_elements(flux)    

    if median(flux) le 0. then fluxbin = 5.e-15 $
    else fluxbin = median(flux)/div
    while (max(flux,iimx,min=mn,subscript_min=iimn)-mn)/fluxbin gt sztol $
      do begin
        if abs(flux[iimx]-fluxbin) gt abs(mn-fluxbin) then begin
            ;;eliminate max flux as larger discrepency
            gd = where(flux ne flux[iimx],ngd)
            flux = flux[gd]
            error = error[gd]
            wave = wave[gd]
        endif else begin
            ;;eliminate min flux as larger discrepency
            gd = where(flux ne flux[iimn],ngd)
            flux = flux[gd]
            error = error[gd]
            wave = wave[gd]
        endelse 
        if not keyword_set(silent) then begin
           if size(spec,/type) eq 7 then $
              print,'measureSN: flux range too large; trim flux and error ',$
                    spec $
           else print,'measureSN: flux range too large; trim flux and error'
        endif                                         ;/silent
    endwhile 

    if median(error) le 0 then errorbin = 5.e-16 $
    else errorbin = median(error)/div
    while (max(error,iimx,min=mn,subscript_min=iimn)-mn)/errorbin gt sztol $
      do begin
        if abs(error[iimx]-errorbin) gt abs(mn-errorbin) then begin
            ;;eliminate max error as larger discrepency
            gd = where(error ne error[iimx],ngd)
            flux = flux[gd]
            error = error[gd]
            wave = wave[gd]
        endif else begin
            ;;eliminate min error as larger discrepency
            gd = where(error ne error[iimn],ngd)
            flux = flux[gd]
            error = error[gd]
            wave = wave[gd]
        endelse 
        if not keyword_set(silent) then begin
           if size(spec,/type) eq 7 then $ 
              print,'measureSN: error range too large; trim flux and error ',$
                    spec $
           else print,'measureSN: error range too large; trim flux and error'
        endif                   ; /silent
    endwhile 

;;Limit to positive flux and corresponding error
;;Flux <= 0. is non-physical
    gd = where(flux gt 0.,ngd)  ;ge 0. causes funniness
    if ngd gt 0 then begin
       ;;center the histograms
       nbins = ceil((max(flux[gd],min=mn)-mn)/fluxbin)
       mn = median(flux[gd]) - 0.5*nbins*fluxbin
        fluxhist = histogram(flux[gd],binsize=fluxbin,locations=fluxloc,min=mn)
        if n_elements(fluxloc) le 3 then begin
            ;; Arrays must be big enough to fit
            fluxhist = [0,0,0,fluxhist,0,0,0]
            fluxloc = (lindgen(n_elements(fluxhist))-3)*fluxbin + fluxloc[0]
            if not keyword_set(silent) then $
              print,'measureSN: buffered flux histogram'
        endif 

       nbins = ceil((max(error[gd],min=mn)-mn)/errorbin)
       mn = median(error[gd]) - 0.5*nbins*errorbin
        errorhist = histogram(error[gd],binsize=errorbin,$
                              locations=errorloc,min=mn)
        if n_elements(errorloc) le 3 then begin
            ;; Arrays must be big enough to fit
            errorhist = [0,0,0,errorhist,0,0,0]
            errorloc = (lindgen(n_elements(errorhist))-3)*errorbin + errorloc[0]
            if not keyword_set(silent) then $
              print,'measureSN: buffered error histogram'
        endif 
        ;; Median values and deviation will be from fitted Gaussian
        ;; (b/c of measure_errors weighting, mdflux/error exactly the
        ;; same) 
        !quiet = 1 ;suppress messages
        gfit = gaussfit(fluxloc,fluxhist,coeff,nterms=3,$
                        measure_errors=sqrt(fluxhist))
        mdflux = coeff[1]
        devflux = coeff[2]
        gfit = gaussfit(errorloc,errorhist,coeff,nterms=3,$
                        measure_errors=sqrt(errorhist))
        mderror = coeff[1]
        deverror = coeff[2]
        snr = mdflux/mderror
        devsnr = snr*sqrt((devflux/mdflux)^2+(deverror/mderror)^2)
    endif else begin
        if not keyword_set(silent) then begin
           if size(spec,/type) eq 7 then $
              print,'measureSN: flux <= 0. for ',spec $
           else print,'measureSN: flux <= 0.'
        endif                   ; /silent
        mdflux = -999.99
        mderror = -999.99
        snr = -999.99
        devflux = -999.99
        deverror = -999.99
        devsnr = -999.99
    endelse
;stop
    if keyword_set(plot) then $
      x_splot,fluxloc,fluxhist,xtwo=errorloc,ytwo=errorhist,/block,$
              psym1=10,psym2=10,xmnx=[min([errorloc,fluxloc],max=mx),mx]
    delvarx,fluxhist,errorhist


endif else begin
;;;;;;;;;;;;;;;
;; Continuum fit and eliminate absorption features
;;;;;;;;;;;;;;;
    safe = {wave:wave,flux:flux,error:error}

    autocontfit,wave,flux,error,buffer=2,binsize=7,region=region,$
                /refin,/silent,view=plot

    mdflux = median(flux[region],/even)
    mderror = median(error[region],/even)
    snr = mdflux/mderror

    devflux = stddev(flux[region])
    deverror = stddev(error[region])
    devsnr = snr*sqrt((devflux/mdflux)^2+(deverror/mderror)^2)

endelse                         ;/cont 

!QUIET = 0                      ;continue suppressing errors

;; Restore arrays
if inflg eq 4 then begin
   wave = svspec.wave
   spec = svspec.flux
   err = svspec.error
endif 

dev = [devflux,deverror,devsnr]
return,[mdflux,mderror,snr]

end
