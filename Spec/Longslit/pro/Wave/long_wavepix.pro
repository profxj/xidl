;+
; NAME:
;   long_wavepix
;
; PURPOSE:
;   Find and trace arc lines
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;  ARCTRC_POS -- Allows for a gap at center (B&C at Dupont!)
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
;------------------------------------------------------------------------------
function long_wavepix, arcimg, tset_slits $
                       , nsig = nsig, box_radius = box_radius $
                       , med_err = med_err1, slit_margin = slit_margin $
                       , VERBOSE = VERBOSE, PKWDTH = PKWDTH, TOLER = TOLER $
                       , FWHM = FWHM1, CHK = CHK, PIXIMG_IN = PIXIMG_IN $
                       , SIG_THRESH = SIG_THRESH, ISLIT = ONLY_SLITS1 $
                       , NOT_THIN = NOT_THIN, arc_ncoeff = arc_ncoeff1 $
                       , ARCTRC_POS=arctrc_pos, FWEIGHT = FWEIGHT1
if N_PARAMS() LT 2 then begin
     print, 'Syntax: wave_pix = long_wavepix(arcimage, left_edge, right_edge, [thresh=, radius=, med_err= ]'
      return, 0
   endif

   ;----------
   ; Set defaults

   IF NOT KEYWORD_SET(NSIG) THEN NSIG = 30.0D
   IF NOT KEYWORD_SET(NAVE) THEN nave = 5L
   if (NOT keyword_set(slit_margin)) then slit_margin = -1
   IF NOT KEYWORD_SET(FORDR) THEN FORDR = 9L
   IF n_elements(FWEIGHT1) GT 0 THEN FWEIGHT = FWEIGHT1 ELSE FWEIGHT = 0
;;  Default values for TRACE_CRUDE
   if (NOT keyword_set(nave)) then nave = 5
   if (NOT keyword_set(box_radius)) then box_radius = 5.
   if (NOT keyword_set(maxshifte)) then maxshifte = 3.0
   if (NOT keyword_set(maxshift0)) then maxshift0 = 3.0
   if (NOT keyword_set(maxerr)) then maxerr = 1.0
   IF NOT KEYWORD_SET(NOT_THIN) THEN THIN = 1
   IF KEYWORD_SET(arc_ncoeff1) THEN arc_ncoeff = arc_ncoeff1 $
   ELSE arc_ncoeff = 3


   ;; ------
   ;; Expand slit set to get left and right edge
   traceset2xy, tset_slits[0], rows, left_edge
   traceset2xy, tset_slits[1], rows, right_edge
   dims = size(arcimg, /dimens)
   nx = dims[0]
   ny = dims[1]
   
   edge_sep = right_edge - left_edge

   ;----------
   ; Extract the arc spectrum using boxcar extraction.
   ; We extract both the left and right half of the slit

;   arc_spec = extract_asymbox2(arcimg, left_edge, right_edge)
;   arc_left = extract_asymbox2(arcimg, left_edge, 0.5*(left_edge+right_edge))
;   arc_right = arc_spec - arc_left

   ; Compute the median slit width for each slit, and analyze the
   ; widest slit first.
   med_width = djs_median(edge_sep,1)
   slit_order = reverse(sort(med_width))
   nslit = n_elements(slit_order)

   if KEYWORD_SET(med_err1) THEN BEGIN
      IF n_elements(med_err1) EQ 1 THEN med_err = replicate(med_err1, nslit) $
      ELSE IF n_elements(med_err1) EQ nslit THEN med_err = med_err1 $
      ELSE message, 'med_err must either be a scale or vector with size nslit'
   ENDIF ELSE med_err = replicate(0.1, nslit)
   
   IF KEYWORD_SET(ONLY_SLITS1) THEN ONLY_SLITS = ONLY_SLITS1 $
   ELSE ONLY_SLITS = lindgen(nslit) + 1L
   IF n_elements(fwhm1) EQ 0 THEN FWHM = replicate(5.0, nslit) $
   ELSE IF n_elements(fwhm1) EQ 1 THEN FWHM = replicate(fwhm1, nslit) $
   ELSE IF n_elements(fwhm1) EQ nslit THEN fwhm = fwhm1 $
   ELSE message, 'Bad size for fwhm input'

   IF KEYWORD_SET(VERBOSE) THEN $
     splog, 'Doing slits in decreasing width order: ', slit_order + 1

   IF KEYWORD_SET(SIG_THRESH1) THEN BEGIN
       IF n_elements(sig_thresh1) EQ 1 THEN $
         sig_thresh = replicate(sig_thresh1, nslit) $
       ELSE IF n_elements(sig_thresh1) EQ nslit THEN sig_thresh = sig_thresh1 $
       ELSE message, 'ERROR: sig_thresh must have either 1 or nslits elements'
   ENDIF ELSE sig_thresh = replicate(3.0, nslit)
     
   tset2d = replicate( {  $
        func    : 'legendre', $
        slit_id :  -1L, $
        dims    : dims, $
        xmin    : 0.0, $
        xmax    : ny-1.0, $
        ymin    : 0.0, $
        ymax    : 1.0, $
        coeff2d : fltarr(5, 5),  $
        coeff2d_inv : fltarr(5, 5) }, nslit)

     ; 5 is the current maximum order of the 2-d fit

   for ii = 0L, nslit-1L do begin 
       jslit = slit_order[ii]
       ionly = WHERE(ONLY_SLITS EQ (jslit+1L), nonly)
       IF nonly EQ 0 THEN CONTINUE
;       jslit = 7L ;; DEBUG
       IF KEYWORD_SET(VERBOSE) THEN $
         splog, 'Working on slit #', jslit+1, ' (', ii+1, ' of ', nslit, ')'
       
      ; The following selects this slit from the image and rectifies
      ; it to a rectangular image.  It works from the slit boundaries
      ; +/- SLIT_MARGIN extra pixels.
      nr = 2*long(ceil(med_width[jslit])/2) 
      ns = (nr + 2*slit_margin + 1) > 1
      ss = (findgen(ns)-slit_margin)/nr 
      yrow = findgen(ny) # (ss*0+1)
      lcen  = ss ## edge_sep[*,jslit] + left_edge[*,jslit] # (ss*0+1)

      arc_rect = interpolate(arcimg, lcen, yrow, cubic=-0.5)

      ;; Approximate the errors...should do this correctly!???
      arc_rect_ivar = 1./(abs(arc_rect + 100) > 1) * (arc_rect GT 1)
      thresh1 = median(arc_rect)*sig_thresh[jslit]
      ;; Identify + trace the arc lines, and fit by polynomials
      arc_rect = reform(arc_rect, ny, ns) ; case where ns=1
      arc_rect_ivar = reform(arc_rect_ivar, ny, ns) ; case where ns=1
      arc1d = fltarr(ny)

      if not keyword_set(ARCTRC_POS) then arctrc_pos = 0.5
      mid = round(ns*arctrc_pos)
      ;; Median filtering is more robust against cosmics
      FOR j = 0L, ny-1L DO BEGIN
          left  = floor(mid - box_radius)
          right = ceil(mid +  box_radius)
          sub_arc = arc_rect[j, left:right]
          djs_iterstat, sub_arc, median = median, sigrej = 2.0
          arc1d[j] = median
       ENDFOR
      IF NOT KEYWORD_SET(PIXIMG_IN) THEN BEGIN
          ;; If we dont have a previous piximg (i.e. from an arc spectrum)
          ;; then we will determine the average trace of the lines first using
          ;; trace_crude, then average this to determine the average profile
          ;; then use that average profile as a crutch for tracing, and the 
          ;; tracing and determination of the average is iterated twice. 
;;         chk = 1
         colors = getcolor(/load)
         x_fndpeaks, smooth(arc1d, 3), peak_pos0, NSIG = NSIG $
                     , PKWDTH = pkwdth, TOLER = TOLER $
                     , THIN = THIN, NORDB = fordr, fweight= fweight
         
     ;;;;;;;;;;; If needed, drop the threshold
         sigvec = [25.D, 20.D, 15.D, 10.D, 7.0]
         isig = 0
         while (peak_pos0[0] EQ -1L OR N_ELEMENTS(PEAK_POS0) LE 5 $
                AND isig LE 4) do begin
;            if ii EQ nslit-1L THEN stop
            x_fndpeaks, smooth(arc1d, 3), peak_pos0, NSIG = sigvec[isig], $
                        /silent, PKWDTH = pkwdth, TOLER = TOLER, $
                        THIN = THIN, NORDB = fordr, fweight = fweight
            isig++
         ENDWHILE
         ;; Identify the highest SNR peaks  (JXP -- did someone comment this out??)
         ;;x_fndpeaks, arc1d, peak_pos0, NSIG = 30.0D, /silent $
         ;;           , PKWDTH = pkwdth, TOLER = TOLER $
         ;;            , /THIN, NORDB = fordr, /fweight
         ;;IF peak_pos0[0] EQ -1L OR N_ELEMENTS(PEAK_POS0) LE 4 THEN $
         ;;   x_fndpeaks, arc1d, peak_pos0 $
         ;;               , NSIG = 7.0D, /silent, PKWDTH = pkwdth, TOLER = TOLER $
         ;;               , /THIN, NORDB = fordr, /fweight
         
         if (keyword_set(CHK)) then begin
            print, ii
            max = 1.2*max(djs_median(arc1d, width = 10))
            plot, arc1d, yrange = [-100, max]
            oplot, peak_pos0, (arc1d[peak_pos0] < 0.9*max) $
                   , psym = 1, color = colors.red
         endif
         

         ;; trace_crude them
         peaks = trace_crude(arc_rect, xstart = peak_pos0 $
                              , ystart = mid, yset = yset $
                              , radius = fwhm[jslit], thresh = thresh1 $
                              , nave = nave $
                              , maxshift0 = maxshift0, maxshifte = maxshifte)
          npeaks = n_elements(peak_pos0)
          IF KEYWORD_SET(VERBOSE) THEN $
            splog, 'Found ', npeaks, ' peaks over ', ns, ' rectified rows'
          ;; Now refine the crude traces using flux weighted tracing
          peakfit = peaks
          FOR i = 1L, 5L DO BEGIN
              peaks_recenter = trace_fweight(arc_rect, peakfit, yset $
                                             , radius = fwhm[jslit])
              xy2traceset, yset, peaks_recenter, tset, ncoeff = arc_ncoeff $
                           , yfit = peakfit, funcname = 'poly' $
                           , maxdev = med_err[jslit], /silent
           ENDFOR
          ;; stack the good traces to determine the average trace profile
          ;; which will be used as a crutch for tracing
          median_abs_dev = djs_median(abs(peaks_recenter-peakfit), 1)
          ;med_err = 0.17 ;; JXP KLUDGE  26 Nov 2014 -- REMOVE!
          good_peaks = where(median_abs_dev LT med_err[jslit] AND $
                             (abs(peaks_recenter-peakfit))[ns/2, *] LT $
                             med_err[jslit], ng)
          
          IF ng LE 1 THEN BEGIN
              splog, 'ERROR: not enough lines found ', ng, ' < 1'
              splog, 'Odds are you do not have enough sky lines to trace.'
              splog, 'or this slit is bad'
              splog, 'We will skip this slit'
              tset2d[jslit].slit_id = jslit
              CONTINUE
          ENDIF ELSE IF ng LT 3 THEN BEGIN
              good_peaks = where(median_abs_dev LT 2.0*med_err[jslit] AND $
                                 (abs(peaks_recenter-peakfit))[ns/2, *] LT $
                                 2.0*med_err[jslit], ng)
          ENDIF
          peak_temp = peaks_recenter[*, good_peaks] $
            - replicate(1.0, ns) # peak_pos0[good_peaks]
          peak_avg = djs_avsigclip(peak_temp, 2, sigrej = 2.0)
          xy2traceset, yset[*, 0], peak_avg, peakset, ncoeff = arc_ncoeff $
                       , yfit = tracefit, funcname = 'poly', /silent 
          ;; this is the avg trace shifted to each peak position
          tracecrutch = tracefit # replicate(1.0, n_elements(peak_pos0)) +  $
            replicate(1.0, ns) # peak_pos0
          niter = 12L
          xfit1 = tracecrutch
          ;; Now perform more accurate tracing using the average trace
          ;; as a crutch. 
          FOR i = 1L, niter DO BEGIN
              xpos1 = trace_fweight(arc_rect, xfit1, yset, radius = fwhm[jslit])
              xy2traceset, yset, xpos1, pos_set1, ncoeff = arc_ncoeff, yfit = xfit1 $
                           , maxdev = med_err[jslit], /silent, funcname = 'poly'
          ENDFOR
          xfit2 = xfit1
          FOR i = 1L, niter DO BEGIN
              xpos2 = trace_gweight(arc_rect, xfit2, yset $
                                    , sigma = fwhm[jslit]/2.3548)
              xy2traceset, yset, xpos2, pos_set2, ncoeff = arc_ncoeff, yfit = xfit2 $
                           , maxdev = med_err[jslit], /silent, funcname = 'poly'
          ENDFOR
          ;; iterate the procedure of finding the average crutch with thers
          ;; refined traces.
          median_abs_dev = djs_median(abs(xpos2-xfit2), 1)
          good_peaks = where(median_abs_dev LT med_err[jslit] AND $
                             (abs(xpos2-xfit2))[ns/2, *] LT med_err[jslit], ng)
          IF ng LE 1 THEN BEGIN
             splog, 'ERROR: not enough lines found ', ng, ' < 1'
             splog, 'Odds are you do not have enough sky lines to trace.'
             splog, 'or this slit is bad'
             splog, 'We will skip this slit'
             tset2d[jslit].slit_id = jslit
             CONTINUE
          ENDIF ELSE IF ng LT 3 THEN BEGIN
             good_peaks = where(median_abs_dev LT 2.0*med_err[jslit] AND $
                                (abs(peaks_recenter-peakfit))[ns/2, *] LT $
                                2.0*med_err[jslit], ng)
          ENDIF
          peak_temp = xfit2[*, good_peaks] $
                      - replicate(1.0, ns) # peak_pos0[good_peaks]
          peak_avg = djs_avsigclip(peak_temp, 2, sigrej = 2.0)
          xy2traceset, yset[*, 0], peak_avg, peakset, ncoeff = arc_ncoeff $
                       , yfit = tracefit2, funcname = 'poly', /silent 
          ;; Now identify low sigma peaks for arc line tracing
          ;x_fndpeaks, arc1d, peak_pos1, NSIG = NSIG, /silent $
          ;            , PKWDTH = pkwdth, TOLER = TOLER $
          ;            , /THIN, NORDB = fordr, /fweight
          ;; Shift the new average trace to each peak position and use 
          ;; as crutch for tracing. 
          peak_pos1 = peak_pos0
          tracecrutch2 = tracefit2 # replicate(1.0, n_elements(peak_pos1)) +  $
            replicate(1.0, ns) # peak_pos1
       ENDIF ELSE BEGIN
          ;; Identify peaks for arc line tracing
          ;x_fndpeaks, arc1d, peak_pos1, NSIG = nsig, /silent $
          ;            , PKWDTH = pkwdth, TOLER = TOLER $
          ;            , /THIN, NORDB = fordr, /fweight

          colors = getcolor(/load)
          x_fndpeaks, arc1d, peak_pos0, NSIG = NSIG $
                      , PKWDTH = pkwdth, TOLER = TOLER $
                      , THIN = THIN, NORDB = fordr, fweight = fweight
     ;;;;;;;;;;; If needed, drop the threshold
          sigvec = [25.D, 20.D, 15.D, 10.D, 7.0]
          isig = 0
          while (peak_pos0[0] EQ -1L OR N_ELEMENTS(PEAK_POS0) LE 5 $
                 AND isig LT 4) do begin
             x_fndpeaks, smooth(arc1d, 3), peak_pos0, NSIG = sigvec[isig], $
                         /silent, PKWDTH = pkwdth, TOLER = TOLER, $
                         THIN = THIN, NORDB = fordr, fweight = fweight
             isig++
          ENDWHILE
          
          if (keyword_set(CHK)) then begin
             print, ii
             max = 1.2*max(djs_median(arc1d, width = 10))
             plot, arc1d, yrange = [-100, max]
             oplot, peak_pos0, (arc1d[peak_pos0] < 0.9*max) $
                    , psym = 1, color = colors.red
;               stop
          endif
          peak_pos1 = peak_pos0

          tracecrutch2 = fltarr(ns, n_elements(peak_pos1))
          pix_rect = interpolate(piximg_in, lcen, yrow, cubic = -0.5)
          xarr=replicate(1.0,ns) ## findgen(ny)
          FOR pp = 0L, n_elements(peak_pos1)-1L DO $ 
            FOR is = 0L, ns-1L DO tracecrutch2[is, pp] = $
            interpol(xarr[*, is], pix_rect[*, is], peak_pos1[pp])
      ENDELSE
      niter = 12L
      xfit1 = tracecrutch2
      yset = rebin(findgen(ns), ns, n_elements(peak_pos1))
      FOR i = 1L, niter DO BEGIN
          xpos1 = trace_fweight(arc_rect, xfit1, yset, radius = fwhm[jslit])
          xy2traceset, yset, xpos1, pos_set1, ncoeff = arc_ncoeff, yfit = xfit1 $
                       , maxdev = med_err[jslit], /silent, funcname = 'poly'
      ENDFOR
      xfit2 = xfit1
      FOR i = 1L, niter DO BEGIN
          xpos2 = trace_gweight(arc_rect, xfit2, yset $
                                , sigma = fwhm[jslit]/2.3548)
          xy2traceset, yset, xpos2, pos_set2, ncoeff = arc_ncoeff, yfit = xfit2 $
                       , maxdev = med_err[jslit], /silent, funcname = 'poly'
      ENDFOR
      peaks_recenter = xpos2
      peakfit = xfit2
      median_abs_dev = djs_median(abs(peaks_recenter-peakfit), 1)
      good_peaks = where(median_abs_dev LT med_err[jslit] AND $
                         (abs(peaks_recenter-peakfit))[ns/2, *] LT med_err[jslit], ng)
      ncoeff = (nr LT 20) ? 2 : 3
      wcoeff = ncoeff + 1
      if ng LT wcoeff then begin
          splog, 'ERROR: not enough lines found ', ng, ' <', wcoeff
          splog, 'Odds are you do not have enough sky lines to trace.'
          splog, 'We receommend you rerun long_reduce with SKYTRACE=0'
          splog, 'Proceed at your own risk'
          tset2d[jslit].slit_id = jslit
          continue
      endif
      IF KEYWORD_SET(VERBOSE) THEN $
        splog, 'Keeping ', ng, ' good peaks out of ', n_elements(peak_pos1)
      if ng LT wcoeff then begin
          splog, 'ERROR: not enough good lines found ', ng, ' <', wcoeff 
          splog, '       Relaxing error threshold by factor of 2'
          splog, '       and lowering order of fit'
          good_peaks = where(median_abs_dev LT 2.0*med_err[jslit] AND $
                             (abs(peaks_recenter-peakfit))[ns/2, *] LT $
                             2.0*med_err[jslit], ng)
          ncoeff = ncoeff-1L
          wcoeff = ncoeff+1L
      endif
      
      xerr = median_abs_dev[good_peaks] ## (ss*0+1)

      xmin = 0.
      xmax = (ny-1.)
;      t = 2.0*(peakfit[ns/2,good_peaks] ## replicate(1,ns) - xmin) $
;             /(xmax-xmin) - 1.
      t = 2.0*(peaks_recenter[*,good_peaks] - xmin) $
        /(xmax-xmin) - 1. 

      ymax = 1.0
      ymin = 0.0
      y = 2.0*(ss # replicate(1,ng) - ymin) / (ymax - ymin) - 1.
      
      long_surf_trace2d, t, y, peakfit[ns/2,good_peaks] ## replicate(1,ns), $
        xerr,surffit, nycoeff=ncoeff, ntcoeff=wcoeff, res=res, mask=tset2d_mask

      tset2d[jslit].slit_id = jslit
      tset2d[jslit].xmin    = xmin
      tset2d[jslit].xmax    = xmax
      tset2d[jslit].ymin    = ymin
      tset2d[jslit].ymax    = ymax
      tset2d[jslit].coeff2d[0:wcoeff-1,0:ncoeff-1] = res

      ;---------------------------------------------------------
      ;  Here's an example of a rectified wavelength image
      ;
      waveimg = flegendre(2.0*findgen(ny)/(ny-1.)-1, wcoeff) # $
        (reform(res, wcoeff, ncoeff) # $
         transpose(flegendre(2.0*ss - 1, ncoeff)))
      ;CHK=1
      IF KEYWORD_SET(CHK) THEN BEGIN
          waveqa = arc_rect
          traceqa = fltarr(ns, ng)
          yqa = rebin(findgen(ns), ns, ng)
          xarr = replicate(1.0, ns) ## findgen(ny)
          FOR pp = 0L, ng-1L DO $ 
            FOR is = 0L, ns-1L DO traceqa[is, pp] = $
            interpol(xarr[*, is], waveimg[*, is], peak_pos1[good_peaks[pp]])
          waveqa[traceqa, yqa] = -10000
          xatv, waveqa, min = -20.0, max = 10*djs_median(arc_rect), /block
      ENDIF
      ;----------------------------------------------------------
      ;  Let's save the inverse function for later use in waveimg
      ;

      keep_for_inv = where(tset2d_mask, nk)
      IF KEYWORD_SET(VERBOSE) THEN $
        splog, 'Keeping ', nk, ' of ', n_elements(tset2d_mask), $
        ' centroids for inverse fit'
      
      tinv = (2.0*(peakfit[ns/2,good_peaks] ## replicate(1,ns) - xmin) $
             /(xmax-xmin) - 1.)[keep_for_inv]
      long_surf_trace2d, tinv, y[keep_for_inv], $
            (peaks_recenter[*,good_peaks])[keep_for_inv], $
            xerr[keep_for_inv], surfinv, $
            nycoeff=ncoeff, ntcoeff=wcoeff, res=resinv

      tset2d[jslit].coeff2d_inv[0:wcoeff-1,0:ncoeff-1] = resinv

      ;---------------------------------------------------------

      ; test chosen radius value 

;      full_radius = djs_median(extract_boxcar(arc_rect, peakfit[*, good_peaks] $
;                                              , radius = box_radius), 1)
;      half_radius = djs_median(extract_boxcar(arc_rect, surffit[*, good_peaks] $
;                                              , radius = box_radius/2.), 1)
;      double_radius = djs_median(extract_boxcar(arc_rect $
;                                                , peakfit[*, good_peaks] $
;                                                , radius = box_radius*2.), 1)
      
;      IF KEYWORD_SET(VERBOSE) THEN BEGIN
;          splog, 'Radius values:', median(half_radius/full_radius), $
;                 median(full_radius/double_radius) 
;          if median(full_radius/double_radius) LT 0.80 then $
;            splog, 'WARNING: Radius: ', radius, ' may be set too small' ;
;          
;          if median(half_radius/full_radius) GT 0.85 then $
;            splog, 'WARNING: Radius: ', radius, ' may be set too large'
;      ENDIF
  endfor



   return, tset2d

end
;------------------------------------------------------------------------------
