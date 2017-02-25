;+ 
; NAME:
; x_fndpeaks   
;    Version 1.1
;
; PURPOSE:
;    Finds a set of peaks in a spectrum.  Generally used with arc line
;    spectra, but valueable in many other areas.
;
; CALLING SEQUENCE:
;   x_fndpeaks, spec, center, NSIG=, PEAK=, /ALL, /SILENT, $
;               PKWDTH=, /THIN, /NORECENT, /FORCE, $
;               FRACPK=, EDGES=, NORDB=, NEDG=, AFUNC=,$
;               MSK=, AUTOFIT=, ICLSE=
;
; INPUTS:
;   spec       - Input spectrum (1D)
;
; RETURNS:
;   center     - Array of peak centers
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  PKWDTH   - Width of peak to center on (pixels) [default: 3L]
;  NORDB    - Order for b-spline fit to the continuum of the
;                spectrum [default: 31L]
;  NSIG=    - Number of sigma signifcance for the peak [default: 5.]
;  NEDG=    - Number of pixels peak must be from the edge of the
;             spectrum [default: 5L]
;  AFUNC=   - Name of function to fit to continuum [default: 'BSPLIN']
;  MSK=     - Mask of good peaks
;  ICLSE=   - Minimum separation between peaks [default: 4L]
;  /ALL     - Pass back all peaks, not just the good ones
;  /THIN    - Allow peaks to be 3 pixel thin [default: 5 pixel]
;  /FORCE   - Force centering in x_centspln
;  /NORECNT - Do not bother to try to recenter the peaks
;  FRACPK=  - Fraction of peak to calculate centroid [default: 0.33]
;  /FWEIGHT - Centroid using flux weighting
;   /FGAUSS   - Centroid arc lines using a gaussian fitting scheme
;
; OPTIONAL OUTPUTS:
;   PEAK=   - Integer values of centeroids of the peaks
;  EDGES=   - Positions of the edges of each peak (usually 0.33 of
;
; COMMENTS:
;
; EXAMPLES:
;   x_fndpeaks, spec, center
;
;
; PROCEDURES/FUNCTIONS CALLED:
; x_centspln
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fndpeaks, spec, center, NSIG=nsig, PEAK=peak, ALL=all, SILENT=silent, $
                PKWDTH=pkwdth, THIN=thin, NORECENT=norecent, FORCE=force, $
                FRACPK=fracpk, EDGES=edges, NORDB=nordb, NEDG=nedg, $
                AFUNC=afunc, FWEIGHT=fweight, TOLER = TOLER, $
                MSK=msk, AUTOFIT=autofit, rms = rms, ICLSE=iclse, $
                FGAUSS=fgauss


;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_fndpeaks, spec, center, NSIG=, /THIN, /NORECENT, NORDB=, PKWDTH=, NEDG= '
    print, '   MSK=, AFUNC=, ICLSE=, /AUTOFIT, /ALL, /FORCE, FRACPK=, EDGES=, PEAK='
    print, '   /FWEIGHT [v1.1]'
    return
  endif 

  npix = n_elements(spec)

; Optional Keywords

  if not keyword_set( TOLER ) then toler = 1.
  ;if not keyword_set( NORDB ) then nordb = 31L  -- JXP 23 Apr 2014
  if not keyword_set( NORDB ) then nordb = 5L
  if not keyword_set( NSIG ) then nsig=5.
  if not keyword_set( PKWDTH ) then pkwdth=3L
  if not keyword_set( NEDG ) then nedg = 5
  if not keyword_set( AFUNC ) then afunc = 'BSPLIN'
  if not keyword_set( MSK ) then msk = bytarr(npix)+1B
  if not keyword_set( ICLSE ) then ICLSE = 4L
  xdat = findgen(npix)
  ; Fit the spectrum with a crude BSPLIN
  if n_elements(AUTOFIT) EQ 0 then begin
     ;; First use this B-spline hack below to create a first guess at
     ;; smooth continuum 

     ;; Modified on August 6, 2013 following JFH modification to
     ;; x_identify
     bkspace = 50.0
     if not keyword_set(NOAMED) then begin ;; JXP Modification to 'help'
        rmed = median(spec, round(bkspace))
                                ;x_splot, arc1d[*,iorder], ytwo=rmed, /bloc
        gdr = where(rmed NE 0.)
                                ;x_splot, arc1d[gdr,iorder]/rmed[gdr], /bloc
        djs_iterstat, spec/rmed[gdr], median = dmed, sigma = dsig 
        gda = where(spec LT rmed+3*dsig)
                                ;x_splot, arc1d[*,iorder], xtwo=gda, ytwo=arc1d[gda,iorder], /bloc
        flg = lonarr(n_elements(spec))
        flg[gda] = 1
     endif else flg = replicate(1, n_elements(xdat))
     spec_set = bspline_iterfit(xdat, spec $
                                , invvar = (spec LE 1d6 and flg and abs(spec) GT 1e-3 )$
                                , bkspace = bkspace $
                                , yfit = autofit_guess $
                                , upper = upper, lower = lower, nord = 3 $
                                , maxrej = 10, outmask = outmask $
                                , /silent, /sticky)
     ;; make this liberal peak finding for the masking
     x_fndpeaks, smooth(spec, 3), peak_pos0, NSIG = 5.0 $
                 , /silent, PKWDTH = pkwdth, TOLER = TOLER $
                 , THIN = THIN, NORDB = NORDB, fweight = fweight, autofit = autofit_guess
     mask = lonarr(npix) + 1L
     IF n_elements(peak_pos0) GT 0 THEN mask[peak_pos0] = 0
     mask_sm = smooth(mask, pkwdth*4.0) ;; Make a liberal mask around each peak 
     ifit = WHERE(mask_sm, nfit)
     IF nfit EQ 0 THEN ifit = lindgen(npix) ;; if there is nothing to fit, fit it all????
     xfit = xdat[ifit] 
     yfit = spec[ifit]
     coeff = robust_poly_fit(xfit, yfit, nordb, /double) 
     autofit = poly(xdat, coeff) 
     ;gd = where(spec NE 0., ngd)
     ; tmpfit = x_fitrej(xdat[gd], spec[gd], afunc, nordb, MSK=msk[gd], $
     ;                   ivar=replicate(1., ngd), $
     ;                   hsigma=2., lsigma=4., rms=rms, FITSTR=fitstr, $
     ;                   /SILENT)
      ;autofit = x_calcfit(xdat, FITSTR=fitstr)
      ;if ptr_valid(fitstr.ffit) EQ 1 then ptr_free, fitstr.ffit
  endif

  ;else begin
  ;    ;; Calculate rms
  ;    if not keyword_set(RMS) then $
  ;      djs_iterstat, spec-autofit, sigma=rms, sigrej=2.0
  ;endelse
  if not keyword_set(RMS) then begin
     gdp = where(abs(spec) GT 1e-3) ;; Another JXP kludge on 2012 Aug 7
     djs_iterstat, (spec-autofit)[gdp], sigma = rms, sigrej = 2.0
  endif
  ; Find all nsig sigma peaks and avoid edges
  gdpix = where( spec GT autofit+nsig*rms AND $
                 xdat GE nedg AND MSK EQ 1B AND $
                 xdat LE npix-nedg-1, ngpix)
  ; Error catch
  if ngpix LE 0 then begin
      if not keyword_set( SILENT ) then $
        print, 'x_fndpeaks: No', nsig, ' sig features!'
      center = -1
      return
  endif

  ; Include only the inflection points (min of 5 pts)
  npk = 0L
  peak = lonarr(npix)
  if not keyword_set( THIN ) then begin
      for i=0L,ngpix-1 do begin
          if ( spec[gdpix[i]-1] GE spec[gdpix[i]-2] AND $
               spec[gdpix[i]] GE spec[gdpix[i]-1] AND $
               spec[gdpix[i]] GE spec[gdpix[i]+1] AND $
               spec[gdpix[i]+1] GE spec[gdpix[i]+2] ) then begin
              peak[npk] = gdpix[i]
              npk = npk + 1
          endif
      endfor
  endif else begin
      for i=0L,ngpix-1 do begin
          if ( spec[gdpix[i]] GE spec[gdpix[i]-1] AND $
               spec[gdpix[i]] GE spec[gdpix[i]+1] ) then begin
              peak[npk] = gdpix[i]
              npk = npk + 1
          endif
      endfor
      ;; Contract all peaks within 2 of one another
      pkmsk = bytarr(npk) 
      spk = peak[sort(peak[0:npk-1])]
      shft2 = shift(spk,-1)
      max = abs(shft2 - spk)
      for i=0L,npk-2 do begin
          if max[i] LE ICLSE then begin
              spk[i] = (spk[i]+spk[i+1])/2.
              pkmsk[i+1] = 1B
          endif
      endfor
      
      ;; Reset
      gd = where(pkmsk EQ 0B, npk)
      peak = spk[gd]
   endelse

  ; Error catch
  if npk LE 0 then begin
      if not keyword_set(silent) then print, 'x_fndpeaks: No peaks found!'
      center = -1
      return
  endif

  peak = temporary(peak[0:npk-1])

  dume = fltarr(2)
  ; EDGES
  if arg_present(EDGES) then edges = fltarr(npk,2)

  if keyword_set( NORECENT ) then begin
      center = double(peak)
      return
  endif else begin
      ; Center up the Peaks 
      center = dblarr(npk)
      tst = dblarr(npk)
      mskc = bytarr(npk)
      for i=0L,npk-1 do begin
          dume[*] = 0.
          ; Worry about limits
          pmn = (peak[i]-pkwdth) > 0
          pmx = (peak[i]+pkwdth) < (npix-1)
          ; Center up on the peak
          if keyword_set( FWEIGHT ) then begin
             ;; JFH 03-27-2014 This code below for flux weighted
             ;; centroiding is majorly flawed. It always defaulted to
             ;; a 3 pixel radius for the centroiding, as the radius
             ;; below was unspecified. Note in cases for which the
             ;; pkwdth was exactly equal to radius, the trace_fweight
             ;; code always returned the midpoint of the input x-values, thus
             ;; ruining the centroiding. I fixed this bug below with
             ;; code which now flux centroids over the range pkwdth
             
             ;;; The following assumes xdat is regularly spaced
             ;;center[i] = x_centfwgt(xdat[pmn:pmx], $
             ;;                       spec[pmn:pmx]-autofit[pmn:pmx], $
             ;;                       radius, SILENT = silent)

             ;; JFH 03-27-2014. My fix below to the bug introduced by
             ;; the x_centfwght routine
             xcen = peak[i]
             FOR k = 0, 19 DO $
                xcen = trace_fweight(spec, xcen, 0L, radius = pkwdth, xerr = xerr $
                                     , invvar = ivar)
             center[i] = xcen
              tst[i] = x_centspln(xdat[pmn:pmx], $
                                     spec[pmn:pmx]-autofit[pmn:pmx], $
                                     fracpk, SILENT=silent, FORCE=force, $
                                     EDGES=dume)
          endif else if keyword_set(FGAUSS) then begin
              guess = [spec[peak[i]]-autofit[peak[i]], peak[i], pkwdth/3., 0., 0.]
              fit = x_gaussfit(xdat[pmn:pmx], $
                               spec[pmn:pmx]-autofit[pmn:pmx], values, $
                               ESTIMATES=guess, SIGMA=sig_fit, NTERMS=5, STATUS=status)
;              clr = getcolor(/load)
;              plot, xdat[pmn:pmx], spec[pmn:pmx], psym=10, $
;                    yrange=[0., max(spec[pmn:pmx])*1.2]
;              oplot, xdat[pmn:pmx], autofit[pmn:pmx], color=clr.red
;              oplot, xdat[pmn:pmx], fit+autofit[pmn:pmx], color=clr.blue
;              oplot, replicate(values[1],2), [0., 1e9], color=clr.green, linesty=1
;              wait, 0.3
              if status EQ 0 then center[i] = values[1] else center[i] = -999
;                  center[i] = -999
;                  stop
;              endelse
          endif else begin  ;; Spline centroiding
              ;; The following does not assumes xdat is regularly spaced
              center[i] = x_centspln(xdat[pmn:pmx], $
                                     spec[pmn:pmx]-autofit[pmn:pmx], $
                                     fracpk, SILENT=silent, FORCE=force, $
                                     EDGES=dume)
          endelse
          ; Edges
          if arg_present(edges) then edges[i,*] = dume
          ; Catch bad centroids 
          if abs(center[i]-xdat[peak[i]]) GT TOLER then begin
;              print,  'DIFF = ', abs(center[i]-xdat[peak[i]])
              center[i] = xdat[peak[i]]
              if arg_present(edges) then edges[i,*] = -1.
          endif else mskc[i] = 1B
      endfor
  endelse
;  printcol, center, center-tst
;  stop
  ; Just pass back the good ones
  if not keyword_set( ALL ) then begin
      gdpk = where(center NE -1 AND mskc,ngdpk)
      if ngdpk EQ 0 then begin
          peak = -1
      endif else begin
          center = center[gdpk]
          if keyword_set(edges) then edges = temporary(edges[gdpk,*])
          peak = peak[gdpk]
      endelse
  endif
  return

end
