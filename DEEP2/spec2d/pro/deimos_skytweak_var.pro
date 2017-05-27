;+
; NAME:
;   deimos_skytweak_var
;
; PURPOSE: 
;   Tweaks the wavelength solution of DEIMOS spSlit files by
;   cross-correlating LOCALLY with a template
;   sky emission spectrum, and fits the local shifts to a
;   fifth-order polynomial (decomposed into Legendre polynomials).
;
; CALLING SEQUENCE:
;   new_fit = deimos_skytweak_var(slit,  sset, color, new_wave1d= , $
;                                 old_wave1d= , flux1d= , $
;                                 template= , fitlam= , $
;                                 shift= , chisq=)
;
; INPUTS:
;   slit       - untweaked DEIMOS spSlit structure (wavelength solution
;                must have been computed using the polyflag method).
;   sset       - bspline sset structure corresponding to slit.
;   color      - which side of the chip are we on ('B' or 'R')?
;
; OPTIONAL INPUTS:

;   	
; KEYWORDS:
;
; OUTPUTS:
;   new_fit    - new coefficients for the wavelength solution
;                polynomial fit. (plug this in to slit.lambdax and
;                call lambda_eval(slit) to get the new 2d wavelength
;                solution.) 
; OPTIONAL OUTPUTS:
;   new_wave1d - 1-d tweaked wavelength solution in a regular 0.2A
;                grid.
;   old_wave1d - 1-d untweaked wavelength solution, as above.
;   flux1d     - 1-d flux, corresponding to new_wave1d and,
;                originally, old_wave1d
;   template   - template spectrum corresponding to old_wave1d
;   fitlam     - wavelengths at which the new wavelength solution was
;                fit.
;   shift      - amount by which old wavelength solution was shifted
;                at the fitlam, before refitting.
;
;   serr   - estimated errors in shift
;
;   chisq      - chi-square of new fit for 2-d wavelength array
;
;   errflag    - 1 if routine encountered an error; zero otherwise.
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;
; REVISION HISTORY:
;   Written by BFG Sept. 02
;----------------------------------------------------------------------


function deimos_skytweak_var,  slit,  sset, color, new_wave1d=new_wave, $
               old_wave1d=fullwave, flux1d=fullflux, $
               template=stemp_flux, fitlam=fitlam, shift=shift, $
               serr=shifterr, chisq = chisq2,  errflag=errflag
;restore,  file='../etc/template.sav'
  errflag =  0
  dir = getenv('DEEP_DIR')+'/'
  if dir eq '/' then message, 'you must set $DEEP_DIR'
  restore,  file=dir+'spec2d/etc/template.sav'

;----restores template sky spectrum from HiRes: 
;----wavelength array in temp_wave, flux in temp_flux

;----smooth the HiRes spectrum with a gaussian to make similar to
;----DEIMOS resolution

  sigma=10
  halfwidth=15
  kernel=findgen(2*halfwidth+1)-halfwidth
  kernel=exp(-kernel^2/2/sigma^2)
  kernel=kernel/total(kernel)

  stemp_flux =  convol(temp_flux,  kernel,  /center)

;stemp_flux =  stemp_flux

;print, 'old fit params: ',  slit.lambdax

  sizex = n_elements(slit.flux[*, 0])
  sizey = n_elements(slit.flux[0, *])
  wave2d = lambda_eval(slit,/double)

;---set up 0.2A wavelength grid
  minlambda =  min(sset.fullbkpt, max=max1) > min(temp_wave, max=max2)
  maxlambda =  max1 < max2
  minlambda = minlambda+1
  maxlambda = maxlambda-1       ; kill off spurious end effects
  dlam =  0.2  ; size of shifts in cross correlation
  nlag = 9   ; number of shifts in cross correlation
  npts =  floor((maxlambda-minlambda)/dlam)
  fullwave = findgen(npts)*dlam + minlambda

;---evaluate object bspline on grid
  fullflux =  bspline_valu(fullwave,  sset)
;---test whether bspline is useful:
  if total(fullflux) gt 0. then begin

;---interpolate template onto grid
  iwave = where((temp_wave ge minlambda) and (temp_wave le maxlambda))
  stemp_flux = interpol(stemp_flux[iwave], temp_wave[iwave],  fullwave)



  window =  100.                ; 100A windows
  nregions =  2*(maxlambda-minlambda)/window

  if color eq 'B' then begin 
     sig_thresh =  100.         ; minimum signal required in stemp_wave over window
     peakthresh =  1.5          ; minimum height for a significant peak
  endif
  if color eq 'R' then begin 
     sig_thresh= 400.           ; minimum signal required in stemp_wave over window
     peakthresh = 10.           ;minimum signal required for a significant peak
  endif
  lag =  indgen(nlag) - nlag/2  ;shift values for x-correlation
  wavlag = lag*dlam
  npoly = 3
  cc =  fltarr(nlag, nregions)     ; cross-correlation values for lag
  cc_err =  cc                  ;  estimated error in cross-correlation

  shiftfit =  fltarr(npoly, nregions) ;polynomial coefficients for fit to cc
  fiterr = shiftfit             ; errors on coefficients

  shift =  fltarr(nregions)     ;peak of fit
  shifterr =  shift             ;error in peak finding

  dofit = shift                 ;which regions' shifts have enough signal to trust
  fitlam = shift                ;central wavelength of fit regions
  fitpix = shift                ; central pixel of fit regions
  tmed = shift                  ;median of each window

  for i=0, nregions-1 do begin
;---march across spectrum in 100A windows, overlapping by 50A
     minwave =  minlambda+i*window/2.
     maxwave =  minwave + window
     if maxwave ge maxlambda then begin
        minwave = maxlambda-window
        maxwave = maxlambda
     endif


     iwin = where((fullwave ge minwave) and (fullwave le maxwave))
     fitpix[i] = min(iwin)+(max(iwin) - min(iwin))/2
     wave = fullwave[iwin]
     fitlam[i] = fullwave[fitpix[i]]
     flux = fullflux[iwin]
     tflux = stemp_flux[iwin]
     djs_iterstat,  tflux,  median=tmp
     tmed[i] = tmp
     
;----does the input bspline exist here?
     if total(flux) ne 0. then begin    
;----is there enough signal to do cross-correlation?
        if total(tflux-tmed[i]) ge sig_thresh then begin
           
;----do cross-correlation, fit with a polynomial, find max
           cc[*, i] = c_correlate(flux,  tflux,  lag)
     
           wpeak =  where(tflux gt peakthresh) ;where there's a significant peak
           if wpeak[0] ne -1 then error = 1/(sqrt(total(tflux[wpeak])/2.)) $
           else error = 1.
           cc_err[*, i] = replicate(error,  nlag)

; ----only fit near the xcorr peak
           maxcorr = max(cc[*, i], maxpix)
           mintofit = (maxpix-2) > 0
           if maxpix lt n_elements(wavlag) - 3 then $
             maxtofit = mintofit+4 else begin
              maxtofit = (maxpix+2) < (n_elements(wavlag)-1)
              mintofit = maxtofit-4
           endelse

           ;print, minwave,  maxwave,  cc_err[0]
           shiftfit[*, i] =  poly_fit(wavlag[mintofit:maxtofit], $
                               cc[mintofit:maxtofit, i], npoly-1, $
                               measure_errors=cc_err[mintofit:maxtofit, i], $
                               sigma=err)
           fiterr[*, i] =  err
           shift[i] = - 0.5*shiftfit[1, i]/shiftfit[2, i]
   
           shifterr[i] = abs(shift[i])*sqrt((fiterr[1, i]/shiftfit[1, i])^2 +$
                                            (fiterr[2, i]/shiftfit[2, i])^2)
;---deweight if we're in the noisy part of the template.
           if minwave ge 8950. then shifterr[i] = sqrt(5.)*shifterr[i]        
           
           dofit[i] = 1 
        endif else begin
           shifterr[i] =  1e5
           dofit[i] = 0
        endelse
     endif else begin
        shifterr[i] =  1e5
        dofit[i] = 0
     endelse
  endfor

  whfit =  where(dofit)
  nfit = n_elements(whfit)
;whfit =  whfit[1:nfit-2];drop first and last points.
  fitlam = fitlam[whfit]
  shift = shift[whfit]
  shifterr = shifterr[whfit]
  fitpix = fitpix[whfit]
  cent_row = floor(sizey/2.)
  cent_wave = wave2d[*, cent_row]
  cent_pix = findgen(sizex)

; find interpolated pixel numbers in central row corresponding 
; to fit wavelengths
  fitpix2 = fitlam
  fitpix2 = interpol(cent_pix, cent_wave, fitlam) 


; the new wavelengths for the fit pixels:
  shift_wave = fitlam+shift


  npix = n_elements(fullwave)
  npix2 = sizex

; do new fits for regular grid and for central row of slit.
  xx = fitpix2/(npix2/2.) -1
  degree = 6
  pxx2 =  fitpix2
  wave_fit2 = svdfit(xx, shift_wave, degree,  /double, /legendre, $
                     yfit=pxx2, measure_errors = shifterr, chisq=chisq)
  chisq2 = chisq
  pxx =  fitpix
  wave_fit = svdfit(fitpix/(npix/2.)-1, shift_wave, degree, /double,/legendre, $
                    yfit=pxx, measure_errors = shifterr, chisq=chisq)


; evaluate new 1-d fit
  pixels = indgen(npix)
  new_wave =  polyleg(pixels/(npix/2.) -1 ,  wave_fit)


; subtract off central dlam from 2-d fit zeroth-order term
  wave_fit2[0] =  wave_fit2[0] - slit.dlam[cent_row]
  
  return, wave_fit2
  endif else begin
    print,  'WARNING: no useful bspline found in deimos_skytweak_var.  Wavelength solution not tweaked'
    errflag =  1
    return, slit.lambdax
  endelse
end


