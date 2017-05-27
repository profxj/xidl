;+
; NAME:
;   lin_arcfit_guess
;
; PURPOSE:
;   Determine initial wavelength solution by comparing spectrum to arc spectrum
;
; CALLING SEQUENCE:
;   wset = lin_arcfit_guess( spec, lambda, intensity, color=color, $
;    [ func=func, bestcorr=bestcorr, acoeff=, dcoeff=, nsteps= ] )
;
; INPUTS:
;   spec       - 1-D spectrum
;   lambda     - lambda of arc lines
;   intensity  - Intensity of arc lines
;
; REQUIRED KEYWORDS:
;   color      - 'red' or 'blue'
;
; OPTIONAL KEYWORDS:
;   func       - Name of fitting function; default to 'legendre'
;   acoeff     - central values of coefficents to explore
;   dcoeff     - range (-.5*dcoeff to .5*dcoeff) of values to explore
;   nsteps     - array of steps to use
;
; OUTPUTS:
;   wset       - traceset (pix -> lambda)
;
; OPTIONAL OUTPUTS:
;   bestcorr   - Correlation coefficient with simulated arc spectrum
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; INTERNAL SUPPORT PROCEDURES:
;   tset_struc()
;   lin_arcfit_iter()
;
; PROCEDURES CALLED:
;   traceset2xy()
;   xy2traceset
;
; REVISION HISTORY:
;   18-Nov-1999  Written by D. Schlegel, Princeton.
;                Excised code from FITARCIMAGE.
;   01-Dec-2000  added acoeff, dcoeff, nsteps keywords
;   24-Aug-2001  branched from arcfit_guess, which uses loglambda.
;          This routine assumes lambda is a roughly linear function of pixno.
;
;-
;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
; NSTEPS = total number of steps to explore in the range
;   [COEFF-0.5*DCOEFF, COEFF+0.5*DCOEFF]
;   except for the first coefficient, for which NSTEPS is unused.
; NSMOOTH = Sigma in pixels for Gaussian-smoothing kernal
; DLAG    = Step in lags for cross-correlation


; for chip no. 2
; plot,spec
; oplot,(lam-6370)/.43,i*2,ps=7


function lin_arcfit_iter, spec, lambda, intensity, $
 aset, dcoeff, nsteps, nsmooth=nsmooth, dlag=dlag, bestcorr=bestcorr

   npix = N_elements(spec)
   nline = N_elements(lambda)

   ncoeff = N_elements(aset.coeff)

   nsz = 4*fix(nsmooth) + 1 ; kernal size (odd)
   gausskern = exp( -( ((nsz-1)/2 - findgen(nsz))^2 ) / nsmooth^2 )
   gausskern = gausskern / total(gausskern)

   npad = fix((npix+1)/2) ; Padding on left and right of spectra before
                          ; cross-correlating

   ; Pad, subtract baseline (with median filter), and smooth input spectrum
   ; Also, take the square-root of the intensity.
   speccorr = spec
   medval1 = median(spec[0:50<npix-1]) ; Median value of left of spectrum
   medval2 = median(spec[npix-50>0:npix-1]) ; Median value of right of spectrum
   speccorr = [fltarr(npad)+medval1, speccorr, fltarr(npad)+medval2] ; Pad
   speccorr = speccorr - median(speccorr, 101<npix) > 1 ; Median-filter
   speccorr = convol(speccorr, gausskern, /center, /edge_truncate) ; Smooth
   speccorr = sqrt(speccorr > 1) ; Weight by the square-root of the intensity

   bestcorr = -2.0
   bestcoeff = 0
   coeff_lo = aset.coeff - dcoeff * (nsteps-1)/(2.*nsteps)
   coeff_lo[0] = aset.coeff[0] ; Set this coefficient to its central value,
                               ; e.g. that passed in ASET.  We will be doing
                               ; the cross-correlation left and right of this.
   tempset = aset


   ; Set minimum number of lags to check equal to 10 pixels
   lagmax = dcoeff[0]
   nlag = fix(lagmax/dlag)
   lags = (indgen(nlag) - fix(nlag/2)) * dlag
   splog, 'Searching lags from ', min(lags), ' to ', max(lags), $
    ' spaced ', dlag, ' pixels', format='(a,i5,a,i5,a,i3,a)'

   ; Loop over all coefficients except the first one
   nsteptot = 1
   for ic=1, ncoeff-1 do nsteptot = nsteptot * nsteps[ic]
   for istep=0, nsteptot-1 do begin
      ; Set the coefficients COEFFI for this step number
      ntemp = nsteptot
      itemp = istep
      for ic=1, ncoeff-1 do begin
         ntemp = ntemp / nsteps[ic]
         j = fix(itemp/ntemp)
         tempset.coeff[ic] = coeff_lo[ic] + (dcoeff[ic]/nsteps[ic]) * j
         itemp = itemp - j * ntemp
      endfor

      ; Construct the simulated arc spectrum
      traceset2xy, tempset, xtemp, tlambda   ; test lambda
      model = fltarr(npix+2*npad)

      if (tlambda[1] GT tlambda[0]) then begin ; Ascending wavelengths
         for iline=0, nline-1 do begin
            qless = (lambda[iline] LT tlambda)
            iloc = (where(qless))[0] - 1
            if (iloc GE 1 AND iloc LE npix-2) then begin
               dx = lambda[iline] - tlambda[iloc]
               dpix = tlambda[iloc+1] - tlambda[iloc]
               model[npad+iloc:npad+iloc+1] = model[npad+iloc:npad+iloc+1] $
                + intensity[iline] * [1-dx/dpix, dx/dpix]
            endif
         endfor
      endif else begin ; Descending wavelengths
         for iline=0, nline-1 do begin
            qless = (lambda[iline] GT tlambda)
            iloc = (where(qless))[0] - 1
            if (iloc GE 1 AND iloc LE npix-2) then begin
               dx = tlambda[iloc] - lambda[iline]
               dpix = tlambda[iloc] - tlambda[iloc+1]
               model[npad+iloc:npad+iloc+1] = model[npad+iloc:npad+iloc+1] $
                + intensity[iline] * [1-dx/dpix, dx/dpix]
            endif
         endfor
      endelse

      ; Smooth the model spectrum.
      ; Also, take the square-root of the intensity.
      model = convol(model, gausskern, /center, /edge_truncate)
      model = sqrt(model > 1) - 1 ; Weight by the square-root of the intensity

      ; Test for the cross-correlation being ill-defined...
      if (stddev(speccorr) EQ 0 OR stddev(model) EQ 0) then begin
         corrval = -1
         icorr = 0
      endif else begin
         corrval = max( c_correlate(speccorr, model, lags), icorr)
      endelse

      if (corrval GT bestcorr) then begin
         bestcorr = corrval
         bestlambda = tlambda
         lagbest = lags[icorr]
; DEBUG PLOTS
; splot,speccorr,xr=[npad,npix+npad]
; soplot,shift(model,-lagbest)*mean(speccorr)/mean(model),color='red'
; print,bestcorr,lagbest,tempset.coeff
      endif

      ; Schlegel counter of step number...
      print, format='("Step ",i5," of ",i5,a1,$)', $
       istep, nsteptot, string(13b)

   endfor

   ; Convert to a trace set with XMIN=0, XMAX=NPIX-1
   xy2traceset, dindgen(npix)-lagbest, bestlambda, wset, ncoeff=ncoeff, $
    xmin=0, xmax=npix-1, maxiter=1

   return, wset
end

;------------------------------------------------------------------------------

function lin_arcfit_guess, spec, lambda, intensity,  $
   color=color, func=func, $
   bestcorr=bestcorr, acoeff=acoeff, dcoeff=dcoeff, nsteps=nsteps

   if (NOT keyword_set(func)) then func = 'legendre'

   npix = N_elements(spec)

   ;---------------------------------------------------------------------------
   ; INITIAL WAVELENGTH SOLUTION
   ;---------------------------------------------------------------------------

   ; Give lin_arcfit_iter initial starting point for wavelength solutions.

;   acoeff = [7250, 4096*.43/2, 0, 0] ; chip 2
;   dcoeff = [2000, 4096*.005/2, 2, 0] 
   acoeff = [7100., 4096.*.32/2, 0., 0.] ; chip 1, 1200 l/mm grating
   if color eq 'red' then acoeff[0] = 8400.  ;red side.
   dcoeff = [2000., 4096.*.005/2, 2., 0.] 

   nacoeff = N_elements(acoeff)
   aset = tset_struc(func, nacoeff)
   aset.xmin = 0
   aset.xmax = npix-1
   aset.coeff = acoeff

   ; First search only varies the first 3 coefficients
   splog, 'Searching about coefficients = ', aset.coeff, format='(a,99f9.3)'
   if not keyword_set(nsteps) then nsteps = [1, 20, 1, 1]
   nsteps = [1, 20, 1, 1]
   wset = lin_arcfit_iter(spec, lambda, intensity, $
    aset, dcoeff, nsteps, nsmooth=6.0, dlag=2, bestcorr=bestcorr)

   ; Second search varies 4 coefficients, but narrows search window
   ; for the first two coefficients
   splog, 'Searching about coefficients = ', wset.coeff, format='(a,99f9.3)'
   nsteps = [1, 10, 10, 1]
   dcoeff[0] = 10 ; Narrow search on coefficient #0
   dcoeff[1] = dcoeff[1] / 5. ; Narrow search on coefficient #1
   wset = lin_arcfit_iter(spec, lambda, intensity, $
    wset, dcoeff, nsteps, nsmooth=4.0, dlag=1, bestcorr=bestcorr)

   splog, 'Initial wavelength fit = ', wset.coeff, format='(a,99f9.3)'
   splog, 'Initial correlation = ', bestcorr

   return, wset
end
;------------------------------------------------------------------------------
