;+
;
; NAME
;       peakinfo.pro
;
; PURPOSE
;       The procedure peakinfo takes a 1-dimensional array containing
;       the spatial flux profile for a given slit and determines the
;       location (column) of the object in the slit and estimates a
;       width of the object's profile. The location of the object is
;       determined from the peak of the spatial flux profile. This
;       peak in the profile is determined via multiple techniques
;       (locating the maximum in the profile, centroiding, and fitting
;       a quadratic to the spatial profile). The width (fwhm) is
;       determined by fitting a spline function to the spatial profile
;       and deriving the fwhm from the spline fit. 
;
;
; SYNTAX
;       peakinfo, sprof, [pkrow, fwhm, pk_quad=pk_quad, $
;          pk_cent=pk_cent, profivar=profivar, s2n=s2n, nbuf=nbuf, $
;          window=window, signif=signif, signal=signal]
;
; INPUTS
;       sprof = a slit profile. That is, sprof is a one-dimensional
;               array which contains the flux profile along the
;               spatial direction of a given slit. That is, sprof
;               contains flux values integrated along the
;               corresponding row of the slit (integrated along the
;               spectral direction). sprof is the output produced by the
;               procedure find_object.pro. In find_object.pro, sky
;               lines are masked out allowing the non-masked flux
;               (object flux) along the slit to be integrated in the
;               spectral direction, thereby estimating the spatial
;               profile of the object.
;       pkrow = a variable to be set equal to the peak column in the
;               spatial profile, that is, to be set equal to the
;               location of the object in the slit. This peak row
;               is estimated by simply determining where the maximum
;               in the spatial profile occurs. Note that the edges of
;               the slit are ignored in finding this maximum so that
;               edge effects do not lead to an inaccurate pkrow
;               value. Also note that the pkrow value will be
;               returned in terms of an integer pixel row value. 
;       fwhm = a variable to be set equal to the estimated fwhm of the
;              object's spatial profile in units of columns. This is 
;              determined using the procedure splinefwhm.pro. Note
;              that the fwhm value returned is in terms of integere
;              number of pixel columns. 
;       pk_quad = a variable to  be set equal to the peak column (see
;                 pkrow) as estimated by fitting a quadratic
;                 function to the spatial profile. Using the intial
;                 pkrow value determined from finding the maximum in
;                 the spatial profile, the procedure mcc_polyfit.pro
;                 is employed to fit a second-degree polynomial to
;                 the columns around the peak. From this fit,
;                 pkrow_qdc is derived. Note that pkrow_qdc is a
;                 float, that is, using the quadratic fitting the
;                 peak column (object location in the slit) is
;                 determined to within a fraction of a pixel. 
;       pk_cent = a variable to be set equal to the peak column (see
;                 pkrow) as estimated by centroiding. That is,
;                 pkrow_ctd is estimated from a centroid around the
;                 peak region. 
;       Note that if peakinfo is unable to fit for a reasonable fwhm
;       value, then the pkrow_ctd and pkrow_qdc values will be set to
;       the flagged value -1.
;      
;       profivar = an optional paramater giving the inverse variance
;                  along the spatial profile. If this parameter is
;                  passed, then peakinfo will return multiple pkrow
;                  and fwhm values (a vector containing one value for
;                  each significant peak in the spatial profile). That
;                  is, the output variables will be vectors will
;                  lengths equal to the number of significant
;                  peaks. The significance of the peak is determined
;                  from a signal-to-noise measurement. 
;       s2n = a variable to be set equal to the signal-to-noise values
;             for each local maximum (peak) found in the spatial
;             profile. s2n will only be set if the significance of
;             each peak is estimated (if the user also supplies the
;             profivar parameter.)
;       nbuf = an optional parameter allowing the user to specify the
;              number of columns to exclude at the edge of the
;              slit. The default value is 2.
;       window = an optional parameter allowing the user to specify
;                the pixel range over which to calculate the
;                signal-to-noise of a peak in the spatial profile. The
;                default window size is 5 pixels - meaning 2 pixels on
;                either side of the peak pixel.
;       signif = an optional parameter which allows the user to
;                specify at what signal-to-noise level should a peak
;                be considered significant. This will onyl matter if
;                the user also passes the profivar parameter. The
;                default is a signal-to-noise ratio of 0. 
;       signal = array of total counts in each peak (per 1d pixel in
;                 extracted spectrum - e.g. in counts/hour/pixel)
;       serendip = set this keyword for detecting the spatial
;                  positions of serendips. if serendip >= 1, then the
;                  spatial profile is ivar smoothed with a 3 pixel
;                  window and the continuum is estimated and
;                  subtracted before detecting the peaks.
;
; KEYWORDS
;       None.
;
; OUTPUTS
;       pkrow, fwhm, pk_quad, pk_cent, s2n (see INPUTS)
;
; PROCEDURES CALLED 
;       splinefwhm2.pro
;       mcc_polyfit.pro
;       bsort.pro
;
; COMMENTS
;       Please note the following flagged values. When peakinfo fails
;       to find a peak in the spatial profile it will return the
;       following values: pkrow = -1, fwhm = -1., s2n = 0.,
;       pk_quad=-1., pk_cent=-1.
;
; EXAMPLES
;       None.
;
; HISTORY
;       Created June 15, 2002 by jnewman and mcc. 
;       Revised July 11, 2002 by mcc - improved documentation and
;          added error messages if sprof parameter is not passed.
;       Revised July 23, 2002 by mcc - revised code such that pkrow is
;          taken to not only be the location of the absolute maximum
;          in the spatial profile, but must also be a local max. This
;          provides a more robust manner for excluding edge effects
;          which distort the spatial profiles. Recall that the spatial
;          profile is determined using the routine
;          find_object.pro. Also changed peakinfo so that it will now
;          take all significant peaks and give an array of pkrow,
;          fwhm, pkrow_qtd, and pkrow_ctd values in return. This will
;          happen if the user passes the profivar parameter.
;
;-

pro peakinfo, sprof, pkrow, fwhm, pk_quad=pk_quad, pk_cent=pk_cent, $
              nbuf=nbuf, profivar=profivar, npix=npix, $
              s2n_fwhm=s2n_fwhm, s2n_window=s2n_window, $
              window=window, signif=signif, serendip=serendip, $
              signal=signal

; make sure that the sprof parameter is supplied.
  if n_elements(sprof) eq 0 then begin
      message, 'user must supply spatial profile (sprof)!', /info
      message, 'See calling sequence in DOC_LIBRARY for more info.'
  endif

; check the values of some of the optional parameters.
  if n_elements(window) gt 0 then window = window[0] else window = 5
  if n_elements(serendip) gt 0 then $
    serendip = serendip[0] ge 1 else serendip = 0
  if n_elements(signif) gt 0 then signif = signif[0] else signif = 0.0
; define the number of pixel rows to exclude on each end of the slit
; (in the spatial direction). default value is 2 rows.
  if n_elements(nbuf) gt 0 then nbuf = nbuf[0] else nbuf = 2

; set the output flag value for this routine.
  piflag = -1.0
; initialize an error catching flag.
  badflag = 0

; determine the number of pixel rows in the slit.
  nrows = n_elements(sprof)

; check that if serendip is set, then profivar was passed.
  if n_elements(profivar) ne nrows and serendip then $
    message, 'user must supply profivar for ivarsmoothing the profile!'

; if we are detecting serendips, smooth the spatial profile. but keep
; the old spatial profile because we want to estimate the s/n from it.
  if serendip then begin
      sprof_old = sprof
      sprof = ivarsmooth(sprof, profivar, 3)
  endif

; locate the local maxima in the spatial profile...to be a maximum a
; point must be greater than 2 points on both sides.
  sup1 = shift(sprof, 1)
  sdw1 = shift(sprof, -1)
  sup2 = shift(sprof, 2)
  sdw2 = shift(sprof, -2)
  mxdex = where(sprof gt sup1 and sprof gt sdw1 and $
                sprof gt sup2 and sprof gt sdw2, npk)
; check that a peak was found...if not, go to the error catching
; section of the program.
  if npk eq 0 then begin
      badflag = 1
      goto, zeropk
  endif

; exclude the edge regions from containing a maximum.
  subdex = where(mxdex gt nbuf-1 and mxdex lt nrows-nbuf, npk)
; check that a peak was still found...if not, go to the error catching
; section of the program.
  if npk eq 0 then begin
      badflag = 1
      goto, zeropk
  endif

; take only the peaks away from the slit ends.
  mxdex = mxdex[subdex]

; if the inverse variance values (profivar) for the spatial profile
; were passed by the user, then check the significance of each peak
; and keep only the significant peaks. [otherwise, we will just take
; the highest peak in the profile.]
  if n_elements(profivar) eq nrows and $
    n_elements(npix) gt 0 then begin
; determine the signal-to-noise (s/n) ratio within a window about each
; local maximum. 
      signal = fltarr(npk)
      noise = fltarr(npk)
      ext = (window - 1)/2
      for ii=0,npk-1 do begin
          dex1 = (mxdex[ii] - ext) > nbuf
          dex2 = (mxdex[ii] + ext) < (nrows-nbuf-1)
          if serendip then begin
              cont = mean([sprof[dex1-1],sprof[dex2+1]]) * window
              signal[ii] = total(sprof[dex1:dex2]) - cont
          endif else signal[ii] = total(sprof[dex1:dex2])
          noise[ii] = sqrt(total(1./profivar[dex1:dex2]))
      endfor
      s2n_window = signal / noise
; only take the peaks with s/n >= signif.      
      submxdex = where(s2n_window ge signif, npk)
; if there are no peaks with s/n >= signif, then take the peak with
; the max s/n. 
      if npk eq 0 then begin
          mxval = max(sprof[mxdex], mxsub) 
          pkrow = mxdex[mxsub]
          s2n_window = s2n_window[mxsub]
      endif else begin
; if there are multiple peaks with s/n >= signif, then sort them based
; on s/n.
          mxdex = mxdex[submxdex]
          s2n_window = s2n_window[submxdex]
          sorted = reverse(sort(s2n_window))
          pkrow = mxdex[sorted]
          s2n_window = s2n_window[sorted]
      endelse
  endif else begin
; if profivar and npix were not set by the user, then simply take the
; highest peak in the spatial profile.
      mxval = max(sprof[mxdex], mxsub) 
      pkrow = mxdex[mxsub]
  endelse

; check how many peaks have survived the above code.
  npk = n_elements(pkrow)

; define output vectors.
  fwhm = fltarr(npk) + piflag
  pk_quad = fltarr(npk) + piflag
  pk_cent = fltarr(npk) + piflag

; now loop and determine fwhm, pk_quad, and pk_cent.
  for ii=0,npk-1 do begin
      pkii = pkrow[ii]
; use splinefwhm2.pro to find the fwhm.
      fwhm1 = splinefwhm2( indgen(pkii+1-nbuf), reverse(sprof[nbuf:pkii]) )
      fwhm2 = splinefwhm2( indgen(nrows-pkii-nbuf), sprof[pkii:nrows-nbuf-1] )
; examine the results, fwhm1 and fwhm2. are the values reasonable?
; splinefwhm2 returns 999.0 when the it fails. if one is bad and the
; other is okay, then simply go with the good one. if both are badm
; then set fwhm equal to the flag value of -1. if both are okay, then
; take the average.
      spflag = 999.0
      if fwhm1 eq spflag then fwhm[ii] = fwhm2 $
      else begin
          if fwhm2 eq spflag then fwhm[ii] = fwhm1 $
          else fwhm[ii] = (fwhm1 + fwhm2) / 2.0 
      endelse
      if fwhm1 eq spflag and fwhm2 eq spflag then fwhm[ii] = piflag

; check that the fwhm value is okay. if so, then proceed to determine
; pk_quad and pk_cent.
      if fwhm[ii] gt 0 then begin
; do a quadratic fitting for peak pixel. 
; define a fitting pixel window around the peak.
          rng = pkii < (nrows-1-pkii) < 5
; fit the quadratic using mcc_ployfit.pro.
          peak = sprof[pkii-rng:pkii+rng]
          rows = (findgen(nrows))[pkii-rng:pkii+rng]
          mcc_polyfit, rows, peak, [0, 1, 2], a=a
; check the quality of the fit.
          if finite(a[1]) and finite(a[2]) then $
            pk_quad[ii] = -a[1] / 2.0 / a[2] $ 
          else pk_quad[ii] = piflag

; centroid to find the peak pixel.
          rng = pkii < (nrows-1-pkii) < floor(fwhm[ii]/2.)
          f = sprof[pkii-rng:pkii+rng]
          x = findgen(2*rng+1) - rng
          delx = total(f * x) / total(f)
          pk_cent[ii] = pkii + delx
      endif
    endfor

; finally, go through the list of peaks and clean out any peaks that
; are actually associated with the same general peak/object in the
; spatial profile. use a resolution of 4-8 pixels. and use the
; centroid positions (if they exist) for this comparison.
    pkc = pk_cent
    dex = where(pk_cent eq piflag, dcnt)
    if dcnt gt 0 then pkc[dex] = pkrow[dex]
    pk_sub = 0
    for jj=0,npk-1 do begin
        resolu = fwhm[jj] / 2.0
        resolu = 4 > resolu
        resolu = resolu < 8
        pkdiff = abs(pkc - pkc[jj])
        pkdiff[jj] = 5*nrows
        minval =  min(pkdiff, minsub)
        if (minval lt resolu and minsub gt jj) or $
          (minval ge resolu) then begin
            if jj gt 0 then pk_sub = [pk_sub, jj]
        endif
    endfor

; trim to only the distinct peaks.
    pkrow = pkrow[pk_sub]
    fwhm = fwhm[pk_sub]
    pk_cent = pk_cent[pk_sub]
    pk_quad = pk_quad[pk_sub]
    if n_elements(s2n_window) gt 0 then s2n_window = s2n_window[pk_sub]
    npk = n_elements(pkrow)

; and last but not least, let's determine the s/n per pixel for each
; peak. here use the fwhm as the window size.
    if n_elements(profivar) eq nrows and $
      n_elements(npix) gt 0 then begin
        signal = fltarr(npk)
        noise = fltarr(npk)
        for ii=0,npk-1 do begin
            if fwhm[ii] le 0 then begin
                signal[ii] = 0.0
                noise[ii] = 1.0E30
            endif else begin
                ext = (fwhm[ii])/2
                dex1 = (pkrow[ii]-ext) > nbuf
                dex2 = (pkrow[ii]+ext) < (nrows-nbuf-1)
                if serendip then $
                  signal[ii] = total(sprof_old[dex1:dex2]) $
                else signal[ii] = total(sprof[dex1:dex2])
                noise[ii] = sqrt(total(1.0/profivar[dex1:dex2]))
            endelse
        endfor
        s2n_fwhm = signal / noise / sqrt(npix[pkrow])
    endif

    zeropk: if badflag then begin 
        message, 'No local maximum found in spatial profile!', /info
        pkrow = piflag
        fwhm = piflag
        s2n_window = 0.0
        s2n_fwhm = 0.0
        pk_quad = piflag
        pk_cent = piflag
    endif
    

end







