;+
; NAME:
;   long_flatfield_specillum
;
; PURPOSE:
;  Determine the illumination function from a flat field exposure
;
;  let's consider x to be the spectral dependence with well sampled
;    (2-pixel) variations
;  y should represent slit-position, and the sampling should be sufficient
;    to resolve slit variations
;  poly_terms will be added to test for possible polynomial variations
;    in x and y.
;  xmin, xmax & ymin, ymax should be added to normalize both parameters to 
;    -1 to 1.
;
;  Assume image represents electrons, and bad pixels are set to negative
;   counts.   We need to mask negative values anyway.
;
;  There is still a issue with scattered light, and how to account for
;   it in the slit illumination response.  Assume that scattered light has
;   already been removed.
;
; CALLING SEQUENCE:
;  long_flatfield_specillum, x, y, image, invvar=, $
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   11-Mar-2005  Written by D. Schlegel, LBL
;-  
;------------------------------------------------------------------------------


pro long_flatfield_specillum, x, y, image, spec_set, illum_set $
                              , invvar = invvar, slitwidth = slitwidth $
                              , modelfit = modelfit, ybkpt = ybkpt $
                              , npoly = npoly, TOL_EDG = TOL_EDG, CHK = CHK $
                              , pixfit = pixfit $
                              , DEBUG = IGOOD, FINEBKPT=finebkpt, xsamp = xsamp $
                              , SLITSAMP = SLITSAMP $
                              ,specsmoothsamp = specsmoothsamp, smoothfit=smoothfit
    chk =1
     npix = n_elements(x)
     t0 = systime(1)
     if n_elements(y) NE npix OR n_elements(image) NE npix then begin
       print, 'Syntax: flatfield_specillum, x, y, image,  $ '
       print, '           spec_set, illum_set, [xsamp=, ysamp=]'
       return
    endif
     IF NOT KEYWORD_SET(TOL_EDG) THEN TOL_EDG = 5.0
     ;; This sets smallest scale of illumuniation function features
;     IF NOT KEYWORD_SET(SLITSAMP) THEN SLITSAMP = 40.0D ;; Should scale with binning
     IF NOT KEYWORD_SET(SLITSAMP) THEN SLITSAMP = 5.0D
     IF NOT KEYWORD_SET(SPECSMOOTHSAMP) THEN SPECSMOOTHSAMP = 50.0 
     ;; This sets minimum rms fluctuation (about unity) required
     ;; to compute an illumination function. 
     IF NOT KEYWORD_SET(ILLUM_THRESH) THEN ILLUM_THRESH = 0.03D
;     IF NOT KEYWORD_SET(ILLUM_THRESH) THEN ILLUM_THRESH = 0.10D
     if NOT keyword_set(xsamp) then xsamp = 0.8
;     if NOT keyword_set(ysamp) then ysamp = 0.1
     if NOT keyword_set(invvar) then invvar = 1.0/((image > 1) + 25.0)

     if n_elements(xmin) NE 1 then xmin = min(x)
     if n_elements(xmax) NE 1 then xmax = max(x)
     if n_elements(ymin) NE 1 then ymin = min(y)
     if n_elements(ymax) NE 1 then ymax = max(y)
     clr = getcolor(/load)

     xrange = xmax - xmin
     if (xrange LE 0) then begin
       print, 'No range found in x, giving up'
       return
     endif

     nxbkpt = long(xrange/xsamp) + 1
     xs  = sort(x)
;
;      First pass, perform a fit to the spectrum in log(counts)
;
     
     log_image = alog(image > 1)
     log_ivar = 1.0 * (image GT 1 AND invvar GT 0)
     inmask_log =  (image GT 1 AND invvar GT 0)
     splog, 'Spectral fit of flatfield', npix, ' Pixels to fit'
     logrej = 0.5
     ;; previously was nord=2
     spec_set = bspline_longslit(x[xs], log_image[xs], log_ivar[xs], $
                                 xs*0.+1, everyn = npix/nxbkpt, nord = 4, $
                                 upper = logrej, lower=logrej, maxrej=5, /groupbadpix, $
                                 yfit=specfit1,inmask = inmask_log, $
                                 outmask=outmask_log, /silent)

;     spec1_set = bspline_iterfit(x[xs], log_image[xs], invvar=log_ivar[xs],  $
;                    everyn=npix/nxbkpt, nord=2, $
;                    upper = 0.2, lower=0.2, maxrej=5, /groupbadpix, $
;                    yfit=specfit1, outmask=tempmask, /silent)
     
     specmodel1 = specfit1*0.
     specmodel1[xs] = exp(specfit1) 
     specmask1 = outmask_log*0
     specmask1[xs] = outmask_log

     IF KEYWORD_SET(CHK) THEN BEGIN
        ;; Plot chi residuals versus spectral and slit position
        badpix = WHERE(outmask_log EQ 0 AND inmask_log GT 0, nbad)
        goodpix = WHERE(outmask_log GT 0 AND inmask_log GT 0, ngood)
        yrange = [0.8*min(specfit1), 1.2*max(specfit1)]
        plot, x[xs[goodpix]], log_image[xs[goodpix]], psym = 3 $
              , xrange = [min(x[xs]), max(x[xs])], yrange =  yrange $
              , xstyle = 1, ystyle = 1, xtitle = 'spectral pixel' $
              , ytitle = 'flat field counts'
        oplot, x[xs], specfit1, col = clr.red
        IF nbad GT 0 THEN oplot, x[xs[badpix]], log_image[xs[badpix]], psym = 1 $
                                 , col =clr.green 
        wait, 1.5
     ENDIF
     

;     Now try a slit-illumination fit

     splog, 'Illumination fit of flatfield slit image'
     splog, 'SLITSAMP = ', SLITSAMP
     ys = sort(y)
     npad = 10000
     yresln = SLITSAMP/slitwidth
     ;; median filter by the resolution sampling interval first
     ;; Compute a quick illumination model by down sampling
     filtimg = image[ys]/(specmodel1[ys] > 1)
     isamp = lindgen(npix/10L)*10L
     samp_width = ceil((n_elements(isamp)*yresln))
     illumquick1 = djs_median(filtimg[isamp] $
                             , width = samp_width, boundary = 'reflect')
     statinds = WHERE(y[ys[isamp]] GT 0.1 AND y[ys[isamp]] LT 0.9, nstat)
     moms = moment(illumquick1[statinds])
     mean = moms[0]
     illum_max = max(abs(illumquick1[statinds]/mean-1.0D))
     ;; If the sub-sampled fluctuations are smaller than threshold/3 
     ;; then skip the full median below which is slow. 
     IF illum_max LE illum_thresh/3.0D THEN BEGIN 
         yin = [-0.2D + 0.2*dindgen(npad)/double(npad-1), y[ys] $
                , 1.0D + dindgen(npad)/double(npad-1)*0.2D]
         normin = [replicate(1.0D, npad), replicate(1.0D, npix) $
                   , replicate(1.0D, npad)]
         maskin = lonarr(2*npad + npix) + 1L
         splog, 'illum_max=' + $
                strcompress(string(illum_max, format = '(F7.3)'), /rem)
         splog, 'Subsampled illum fluctuations < threshold/3=' $
                + strcompress(string(100.0*illum_thresh/3.0D $
                                     , format = '(F4.2)'), /rem) + '%'
         splog, 'No illum function applied for this slit'
     ENDIF ELSE BEGIN 
         illumquick = interpol(illumquick1, y[ys[isamp]], y[ys])
         chi_illum = (image[ys]-specmodel1[ys]*illumquick)*sqrt(invvar[ys])
         medmask = (invvar[ys] GT 0.0) AND abs(chi_illum) LT 10.0
         ifit = WHERE(medmask, nmed)
         med_width = ceil((nmed*yresln))
         normimg = djs_median(filtimg[ifit], width = med_width $
                              , bound = 'reflect')
         sig_res = med_width/15L
         nhalf =  long(sig_res)*4L
         xkern = dindgen(2*nhalf+1)-nhalf
         kernel = gauss1(xkern, [0.0, sig_res, 1.0])
         normimg = convol(normimg, kernel, /edge_truncate)
         ;;normimg = interpol(normimg, y[ys[ifit]], y[ys])
         ;;statinds = WHERE(y[ys[ifit]] GT 0.1 AND y[ys] LT 0.9)    
         ;; KHRR CHANGED THE ABOVE LINE TO THE LINE BELOW -- 29 Sep 2014
         statinds = WHERE(y[ys[ifit]] GT 0.1 AND y[ys[ifit]] LT 0.9)
         moms = moment(normimg[statinds])
         mean = moms[0]
         normimg = normimg/mean
         ;; compute median value of normimg edge pixels 
         lmed = djs_median(normimg[0:9])
         rmed = djs_median(normimg[(nmed-10L):(nmed-1L)]) 
         ;; mask regions where illumination function takes on extreme values
         imask = (specmodel1[ys[ifit]] GT 1   AND $
                  invvar[ys[ifit]]     GT 0.0 AND $
                  finite(normimg) EQ 1)
         nanpix = WHERE(finite(normimg) EQ 0, nnan)
         IF nnan GT 0 THEN message, 'Inifinities in normimg'
         illum_max = max(abs(normimg[statinds]-1.0D))
         IF illum_max LE illum_thresh THEN BEGIN 
             yin = [-0.2D + 0.2*dindgen(npad)/double(npad-1), y[ys] $
                    , 1.0D + dindgen(npad)/double(npad-1)*0.2D]
             normin = [replicate(1.0D, npad), replicate(1.0D, npix) $
                       , replicate(1.0D, npad)]
             maskin = lonarr(2*npad + npix) + 1L
             splog, 'illum_max=' + $
                    strcompress(string(illum_max, format = '(F7.3)'), /rem)
             splog, 'Max illum fluctuations < ' + $
                    strcompress(string(100.0*illum_thresh, format = '(F4.2)') $
                                , /rem) + '%'
             splog, 'No illumination function applied for this slit'
         ENDIF ELSE BEGIN 
             yin = [-0.2D + 0.2*dindgen(npad)/double(npad-1), y[ys[ifit]] $
                    , 1.0D + dindgen(npad)/double(npad-1)*0.2D]
             normin = [replicate(lmed, npad), normimg, replicate(rmed, npad)]
             maskin = [replicate(1.0D, npad), double(imask) $
                       , replicate(1.0D, npad)]
         ENDELSE
      ENDELSE
     IF KEYWORD_SET(ybkpt) THEN fullbkpt = ybkpt $
     ELSE begin
        if not keyword_set(FINEBKPT) then bksp = 10. else bksp = 50.
        ybkpt = bspline_bkpts(yin, nord = 4, bkspace = yresln/bksp, /silent)
     endelse
     ;; Fixed a bug here where fits were crashing because of too
     ;; finely space breakpoints 6-8-2013 JFH. Changed yrresln/50 to yresln/10
     ;; JXP -- This broke Kast.  Have added a 'fix'  8-21-2013
     ;;ELSE ybkpt = bspline_bkpts(yin, nord = 4, bkspace = yresln/50.0, /silent)
     illum_set = bspline_longslit(yin, normin, maskin, yin*0.+1.0D $ 
                                  , yfit = illumfit1, upper = 5, lower = 5 $
                                  , outmask = tempmask2 $
                                  , fullbkpt = ybkpt $
                                  , /silent)
     ;;chk = 1
     IF KEYWORD_SET(CHK) THEN BEGIN
         yrange = [0.8*min(illumfit1), 1.2*max(illumfit1)]
         plot, y[ys], image[ys]/(specmodel1[ys] > 1), psym = 3 $
               , xrange = [0.0, 1.0], yrange =  yrange $
               , xstyle = 1, ystyle = 1, xtitle = 'slit fraction' $
               , ytitle = 'median filtered illum func'
         oplot, yin, illumfit1, col = clr.red
         IF KEYWORD_SET(normimg) THEN oplot, y[ys[ifit]], normimg $
           , col = clr.green
         wait, 1.5
     ENDIF
     IF NOT KEYWORD_SET(PIXFIT) THEN RETURN
         
      illummodel1 = y[ys]*0 
      illummodel1[ys] = bspline_valu(y[ys], illum_set) 

      if NOT keyword_set(npoly) then npoly = 7L ;4L  
      profile_basis = fpoly(2.0D*y[xs] - 1.0D, npoly)
      ;poly_basis = fpoly(2.0D*y[xs] - 1.0D, npoly)
      ;IF illum_max LE illum_thresh THEN $
      ;  profile_basis = poly_basis $
      ;ELSE profile_basis = [[(illummodel1[xs])], [poly_basis]]
      splog, 'Illumination+scattered light fit of flatfield slit image'
      ;; pre-mask pixels which are outliers from our current model. 
      ;; As these seem to be breaking the fits
      clip = 15.0d
      inmask = (invvar[xs] GT 0.0) AND $
        abs((image[xs]-specmodel1*illummodel1[xs])*sqrt(invvar[xs])) LT clip
;        (y[xs] GE TOL_EDG/SLITWIDTH)  AND $
;        (y[xs] LE (1.0D - TOL_EDG/SLITWIDTH))

      ;; Create a prelminary normalized flat from the first pass
      ;; separable decomposition = specmodel*illumodel. Now fit out
      ;; more complicated spectrally dependent polynomial
      ;; dependendence using bspline_longslit
      imag_norm_pix = image[xs]/(specmodel1*illummodel1[xs])
      ;; Here we ignore the formal photon counting errors and simply
      ;; assume that a typical error per pixel. This guess is somewhat aribtrary
      ivar_imag_norm = inmask/(0.01^2)
      sigrej_illum = 3.0
      ;; The above presumes that we will mask ~ +-15% outliers from
      ;; this model.
      scatillum_set = bspline_longslit(x[xs], imag_norm_pix $
                                       , ivar_imag_norm $
                                       , profile_basis, yfit = yfit $
                                       ,  bkspace = specsmoothsamp  $
                                       , outmask = outmask $
                                       , inmask = inmask $
                                       , nord = 4 $
                                       , maxrej=10, /groupbadpix $
                                       , upper = sigrej_illum $
                                       , lower = sigrej_illum) ;, DEBUG=[IGOOD,xs])
      IF total(yfit) LE 1.0e-6 THEN yfit = 0.0*xs + 1.0 
;      , upper = 7, lower = 7)
      ;scat2 = bspline_longslit(x, image $
                                       ;, invvar*inmask $
                                       ;, profile_basis, yfit = yfit2 $
                                       ;, fullbkpt = spec_set.fullbkpt $
                                       ;, outmask = outmask $
                                       ;, inmask = inmask $
                                       ;, nord = 4 $
                                       ;, upper = 9, lower = 9)
      
;      IF KEYWORD_SET(DEBUG) THEN BEGIN
;          indx = 0
;          prof_set = $
;            create_bsplineset(scatillum_set.fullbkpt, scatillum_set.nord)
;          prof_set.bkmask = scatillum_set.bkmask
;          prof_set.coeff = scatillum_set.coeff[indx, *]
;          prof_model = bspline_valu(x[xs], prof_set)*profile_basis[*, indx]
;      ENDIF
      ;; Plot chi residuals versus spectral and slit position
      badpix = WHERE(outmask EQ 0 AND inmask GT 0, nbad)
      goodpix = WHERE(outmask GT 0 AND inmask GT 0, ngood)
      residual = (imag_norm_pix - yfit) 
      weight = ivar_imag_norm 
      chi    = residual*sqrt(weight)
      ;CHK=1
      IF KEYWORD_SET(CHK) THEN BEGIN
          plot, x[xs[goodpix]], residual[goodpix], psym = 3 $
                , xrange = [min(x[xs]), max(x[xs])] $
                , yrange = [-0.05, 0.05] $
                , xstyle = 1, xtitle = 'spectral pixel', ytitle = 'residual'
          IF nbad GT 0 THEN oplot, x[xs[badpix]], residual[badpix], psym = 1 $
            , col =clr.green
          wait, 3.0
          plot, y[xs[goodpix]], residual[goodpix], psym = 3, xrange = [0.0, 1.0] $
                , yrange = [-0.05, 0.05], xstyle = 1 $
                , xtitle = 'slit fraction' $
                , ytitle = 'residual'
          IF nbad GT 0 THEN oplot, y[xs[badpix]], residual[badpix], psym = 1 $
                                   , col =clr.green 
          wait, 1.5
       ENDIF
      smoothfit  =  x*0.
      smoothfit[xs] = yfit
      modelfit = x*0
      modelfit[xs] = yfit*specmodel1*illummodel1[xs]
      splog, 'Elapsed time = ', systime(1)-t0, ' sec'
      return
  end
               
                           
