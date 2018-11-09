;+
; NAME:
;   long_objfind
;
; PURPOSE:
;   Find the location of objects within each slit mask
;
; CALLING SEQUENCE:
;   objstruct = long_objfind( image, tset_slits=, $
;    [fwhm=, nperslit=, peakthresh=, skymask=, objmask= ] )
;
; INPUTS:
;   image      - Image for finding the slits, which would typically
;                be a flat-field image, an arc image, or a sum of those
;   tset_slits - Trace sets with slit start/end positions
;
; OPTIONAL INPUTS:
;   fwhm       - FWHM for convolving flux along the slit before peak-finding;
;                default to 3.0 pix.
;                (Also, do not allow peaks closer to each other than FWHM pix.)
;   nperslit   - Do not find more than this many objects per slit; default to 10
;   peakthresh - Flux threshhold for finding objects; the flux must be
;                at least this fraction of the brightest object in each slit;
;                default to 0.
;   absthresh  - Absolute flux threshold for finding objects; the peakflux must 
;                be at least htis large; default to 0. If both peakthresh and 
;                absthresh are set, absthresh overrides peakthresh. 
;   SIG_THRESH - Sigma threshold for objects  [default=5.0]
;   OBJTHRESH  - threshold for object masking
;   /SIMPLE_SUB -- Do the simple sky sub (median) no matter what
;   CLIP_SLIT -- Allow for truncated slits (e.g. ESI order=10)
; OUTPUTS:
;   objstruct  - Structure with object parameters
;
; OPTIONAL OUTPUTS:
;   skymask    - Image of the cleanest sky pixels, =1 for good sky pixels
;
; COMMENTS:
;   Each slit is smashed to a single vector in the spatial direction,
;   taking care to align each wavelength in the spatial direction first.
;   On this smashed flux vector, we then search for significant peaks.
;   We then compute the expected position of each object at each wavelength,
;   by assuming that each object always falls in the same position within
;   the slit.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   extract_asymbox2()
;   find_npeaks()
;   gaussian()
;   long_slits2mask()
;   long_slits2x()
;   splog
;   traceset2xy
;
; INTERNAL SUPPORT ROUTINES:
;   long_obj_create()
;
; REVISION HISTORY:
;   10-Mar-2005  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
; Create the return structure for NOBJ objects
function long_obj_create, nobj, slitid=slitid, ny=ny

   objstr = create_struct( $
    'OBJID' , 0L, $
    'SLITID', 0L, $
    'XFRACPOS', 0.0, $
    'PEAKFLUX', 0.0, $
    'MASKWIDTH', 0.0, $
    'FWHM', 0.0, $
    'FLX_SHFT_WAV', 0.0, $
    'FLX_SHFT_SPA', 0.0, $                   
    'FWHMFIT'  , fltarr(ny), $
    'XPOS'  , fltarr(ny), $
    'YPOS', findgen(ny), $
    'HAND_AP', 0L, $
    'HAND_X', 0.0, $
    'HAND_Y', 0.0, $
    'HAND_MINX', 0.0, $
    'HAND_MAXX', 0.0, $
    'HAND_FWHM', 0.0, $
    'HAND_SUB', 0L)                        
;    'XPOSIVAR'  , fltarr(ny), $
     if (keyword_set(slitid)) then objstr.slitid = slitid
   if (keyword_set(nobj)) then objstr = replicate(objstr, nobj)

   return, objstr
end

;------------------------------------------------------------------------------
; Function to deal with clipped slits (e.g. ESI)
function long_obj_clipslit, img, clip_slit, x, y, objs
  
  xinvvar = x
  xinvvar[*] = 1.
  for yy=0L,n_elements(CLIP_SLIT)-1 do begin
     mt = where(objs.slitid EQ CLIP_SLIT[yy], nslit)
     for ss=0L,nslit-1 do begin
        idx = mt[ss]
        badp = where(abs(img[round(x[*,idx]), round(y[*,idx])]) LT 1e-4, nbadp)
        if nbadp GT 0 then xinvvar[badp,idx] = 0.
     endfor
  endfor
  
   return, xinvvar
end
;------------------------------------------------------------------------------
function long_objfind, image, tset_slits=tset_slits $
                       , fwhm = fwhm1, nperslit = nperslit1 $
                       , ncoeff = ncoeff1 $
                       , peakthresh = peakthresh1, ABSTHRESH = ABSTHRESH1 $
                       , skymask = skymask_new, objmask = objmask $
                       , OBJTHRESH = OBJTHRESH, CRUDE = CRUDE $
                       , sky_fwhm = sky_fwhm, sigma_sky = sigma_sky $
                       , SILENT = SILENT, POS_SET = POS_SET $
                       , PEAK_SMTH=PEAK_SMTH, SIG_THRESH=sig_thresh1 $
                       , HAND_X = HAND_X $
                       , HAND_MINX = HAND_MINX1, HAND_MAXX = HAND_MAXX1 $
                       , HAND_Y = HAND_Y $
                       , HAND_FWHM = HAND_FWHM1, HAND_SUB = HAND_SUB1 $
                       , INVVAR = invvar, SIMPLE_SUB=simple_sub $
                       , VERBOSE = VERBOSE, STDTRACE = STDTRACE, ISLIT = ISLIT $
                       , CLIP_SLIT=clip_slit

  ;;----------
   ; Test inputs

   if (keyword_set(image) EQ 0) then $
    message, 'Must specify IMAGE'

   if NOT keyword_set(invvar) then  invvar = 1.0/(abs(image) + 100.)
   
   ;----------
   ; Set defaults

   dims = size(image, /dimens)
   nx = dims[0]
   if nx LE 2 then begin
     splog, 'Not enough columns in IMAGE', nx
     return, 0
   endif
   ;; Default values for TRACE_CRUDE
   if (NOT keyword_set(ksize)) then ksize = 5
   if (ksize LT 1 OR ksize GE nx/2-1) then $
     message, 'Invalid kernel size KSIZE'
   if (NOT keyword_set(nave)) then nave = 3
   if (NOT keyword_set(maxshifte)) then maxshifte = 0.1
   if (NOT keyword_set(maxshift0)) then maxshift0 = 1.0
   IF NOT KEYWORD_SET(PEAK_SMTH) THEN PEAK_SMTH=5.0
   IF NOT KEYWORD_SET(CRUDE_TOL) THEN CRUDE_TOL = 30.0D
   IF NOT KEYWORD_SET(NCOEFF1) THEN NCOEFF = 5 $
   ELSE NCOEFF = NCOEFF1
   func_crude = tset_slits[0].func
   dim_crude = size(tset_slits.COEFF, /dim)
   ncoeff_crude = dim_crude[0]

   ny = dims[1]
   ymid = ny / 2
   if (keyword_set(fwhm1)) then fwhm = fwhm1 else fwhm = 3.0 ; in pixels
;   if (NOT keyword_set(fwhm)) then fwhm = 3.0 ; in pixels
   if (keyword_set(peakthresh1)) then peakthresh = peakthresh1 $
    else peakthresh = 0.
   if (keyword_set(SIG_THRESH1)) then SIG_THRESH = SIG_THRESH1 $
   else SIG_THRESH=5.0
   if (keyword_set(absthresh1)) then absthresh = absthresh1 $
   else absthresh = 0.
   if (keyword_set(nperslit1)) then nperslit = nperslit1 $
    else nperslit = 10
   IF NOT KEYWORD_SET(OBJTHRESH) THEN OBJTHRESH = 0.5D

   objstruct = 0
   traceset2xy, tset_slits[0], yy1, xx1
   traceset2xy, tset_slits[1], yy2, xx2

   slitmask = long_slits2mask(tset_slits, nslit = nslit)
   skymask = (slitmask GT 0)
   objmask = slitmask*0

   ;; Identify slits for specified aperture positions 
   IF KEYWORD_SET(HAND_X) THEN BEGIN
       IF n_elements(HAND_X) NE n_elements(HAND_Y) THEN $
         message, 'HAND_X and HAND_Y must be set'
;       naps = n_elements(HAND_X)
       handslits = slitmask[HAND_X, HAND_Y]
       IF KEYWORD_SET(HAND_FWHM1) THEN BEGIN
           IF n_elements(HAND_FWHM1) EQ 1 THEN $
             hand_fwhm = replicate(hand_fwhm1, n_elements(HAND_X)) $
           ELSE IF n_elements(HAND_FWHM1) EQ n_elements(HAND_X) THEN $
             hand_fwhm = hand_fwhm1 $
           ELSE IF n_elements(HAND_FWHM1) NE n_elements(HAND_X) THEN $
             message, 'ERROR: Invalid number of elements for HAND_FWHM'
       ENDIF ELSE hand_fwhm = fltarr(n_elements(HAND_X))
       IF KEYWORD_SET(HAND_SUB1) THEN BEGIN
          IF n_elements(HAND_SUB1) EQ 1 THEN $
             hand_sub = replicate(hand_sub1, n_elements(HAND_X)) $
          ELSE IF n_elements(HAND_SUB1) EQ n_elements(HAND_X) THEN $
             hand_sub = hand_sub1 $
          ELSE IF n_elements(HAND_SUB1) NE n_elements(HAND_X) THEN $
             message, 'ERROR: Invalid number of elements for HAND_SUB'
       ENDIF ELSE hand_sub = lonarr(n_elements(hand_x)) ;; If HAND_SUB is not specificed default to not subtract these apertures
       IF KEYWORD_SET(HAND_MINX1) THEN BEGIN
          IF n_elements(HAND_MINX1) EQ 1 THEN $
             hand_minx = replicate(hand_minx1, n_elements(HAND_X)) $
          ELSE IF n_elements(HAND_MINX1) EQ n_elements(HAND_X) THEN $
             hand_minx = hand_minx1 $
          ELSE IF n_elements(HAND_MINX1) NE n_elements(HAND_X) THEN $
             message, 'ERROR: Invalid number of elements for HAND_MINX'
       ENDIF ELSE hand_minx = fltarr(n_elements(HAND_X))
       IF KEYWORD_SET(HAND_MAXX1) THEN BEGIN
          IF n_elements(HAND_MAXX1) EQ 1 THEN $
             hand_maxx = replicate(hand_maxx1, n_elements(HAND_X)) $
          ELSE IF n_elements(HAND_MAXX1) EQ n_elements(HAND_X) THEN $
             hand_maxx = hand_maxx1 $
          ELSE IF n_elements(HAND_MAXX1) NE n_elements(HAND_X) THEN $
             message, 'ERROR: Invalid number of elements for HAND_MAXX'
        ENDIF ELSE hand_maxx = fltarr(n_elements(HAND_X))
   ENDIF
      

   ; Make gaussian kernels
   ngpix = (2 * ceil(1.5*fwhm) + 1) < (nx-1)
   width = 3
   gkern = gaussian(findgen(ngpix) - (ngpix-1)/2, [1.,0.,fwhm/2.305])
   gkern = gkern / total(gkern)
   ;----------
   ; Loop over each slit
   
   IF KEYWORD_SET(ISLIT) THEN BEGIN
       IF islit GT nslit THEN message $
         , 'ERROR: islit not found. islit cannot be larger than nslit'
       nreduce = 1
       slit_vec = [islit]
   ENDIF ELSE BEGIN
       nreduce = nslit
       slit_vec = lindgen(nslit) + 1L
   ENDELSE
   mask = (invvar GT 0.0)
   for jj = 0L, nreduce-1L DO BEGIN
       slitid = slit_vec[jj]
       thisimg = image * (slitmask EQ slitid)*mask
       ximg = long_slits2x(tset_slits, slitid = slitid)

       ;; Smash the image (for this slit) into a single flux vector.
       ;; How many pixels wide is the slit at each Y?
       xsize = xx2[*,slitid-1] - xx1[*,slitid-1]
       ;; How many pixels to sample this at?
       nsamp = ceil(djs_median(xsize))
;      nsamp = ceil(max(xsize)) ; changed by JFH aug 2007. Not sure if will
;;     cause problems but using the maximum nsamp doesn't seem to make sense
       
;     Mask skypixels with 2 FWHM of edge
;
       IF NOT KEYWORD_SET(SKY_FWHM) THEN SKY_FWHM = FWHM
       border = (2*sky_fwhm/nsamp) < 0.03
       skymask = skymask * (1B-((slitmask EQ slitid) AND $
                                ((ximg LT border) OR (ximg GT 1 - border))))
       
       left_asym = rebin(xx1[*,slitid-1], ny, nsamp) + $
                   (xsize/nsamp) # (findgen(nsamp))
       right_asym = left_asym + (xsize/nsamp) # replicate(1,nsamp)
       
       flux_spec = extract_asymbox2(thisimg, left_asym, right_asym)
       fluxvec = djs_avsigclip(flux_spec, 1, sigrej = 4)

;     Changed reject threshold from sigrej=2. To allow for more 
;     dynamic range since some good pixels were being rejected. 
;     Median-filter with a filter much larger than the FWHM
;     don't reflect, it just finds objects near the edge
      if (nsamp LT 9.*fwhm) OR keyword_set(SIMPLE_SUB) then fluxsub = fluxvec - median(fluxvec) $
      ELSE BEGIN 
         fluxsub = fluxvec - djs_median(fluxvec, width = PEAK_SMTH*fwhm)
         ;; This little bit below deals with degenerate cases for which the 
         ;; slit gets brighter toward the edge, i.e. when alignment stars
         ;; saturate and bleed over into other slits. In this case
         ;; the median smoothed profile is the nearly everywhere the 
         ;; same as the profile itself, and fluxsub is full of zeros
         ;; (bad!). If 90% or more of fluxsub is zero, default to use
         ;; the unfiltered case
         isub_bad = WHERE(fluxsub EQ 0.0, nsub_bad)
         frac_bad = double(nsub_bad)/double(nsamp)
         IF frac_bad GT 0.9 THEN fluxsub = fluxvec - median(fluxvec)
      ENDELSE
      ; Convolve this with a gaussian
      fluxconv = convol(fluxsub, gkern, /center, /edge_truncate)
      ; Find the peaks
      xcen = long_find_nminima(-fluxconv, nfind = nperslit $
                               , width = width $
                               , minsep = (fwhm > width), ypeak = ypeak1 $
                               , npeak = npeak)
      ypeak = interpolate(fluxconv, xcen)
      ;; This below masks out the 2.0*FWHM from the stats JFH 09/11
      imask = lonarr(nsamp) + 1L
      xvec = findgen(nsamp)
      FOR zz = 0L, npeak-1L DO BEGIN
         ibad = WHERE(abs(xvec-xcen[zz]) LE FWHM, nbad)
         IF nbad GT 0 THEN imask[ibad] = 0L
      ENDFOR
      ;; Omit edge pixels and pixels within a FWHM of the object 
      ;; from these stats.  JFH 09/11
      igd = WHERE(imask AND xvec GT 3L AND xvec LE (nsamp-3L), ngd)
      IF ngd EQ 0 THEN igd = lindgen(n_elements(fluxconv))
      djs_iterstat, fluxconv[igd], sigma = skythresh, sigrej = 1.5
      djs_iterstat, fluxconv[igd], sigma = sigma, sigrej = 2.5
      IF skythresh EQ 0.0 THEN skythresh = sigma
; jm10may25ucsd; if SKYTHRESH is *still* zero then mask out the pixels
; that are zero
      if skythresh EQ 0.0 THEN begin 
         good = where(fluxconv gt 0.0,ngood)
         if (ngood eq 0L) then message, 'This should not happen'
         djs_iterstat, fluxconv[good], sigma=skythresh, sigrej=1.5
      endif
      
      ; Get rid of peaks within 3% of slit edge
      not_near_edge = where(xcen GT nsamp*0.03 AND xcen LT nsamp*0.97,npeak)
      if npeak GT 0 then begin
        xcen = xcen[not_near_edge]
        ypeak = ypeak[not_near_edge]
      endif
      if npeak GT 1 then begin
         threshold = SIG_THRESH*sigma 
         IF keyword_set(peakthresh) then $
            threshold = threshold > (peakthresh*max(ypeak))
         IF KEYWORD_SET(ABSTHRESH) THEN $
            threshold = threshold > ABSTHRESH
         IF KEYWORD_SET(VERBOSE) THEN $
            splog, 'Slitid: ', slitid, ' Using a threshold of ', $
                   threshold, max(ypeak)
                                ; Trim to only objects above threshhold
         ikeep = where(ypeak GT threshold, npeak)
         if npeak GT 0 then begin
            xcen = xcen[ikeep]
            ypeak = ypeak[ikeep]
         endif
      endif

    IF npeak NE 0 THEN BEGIN
        handaps = lonarr(npeak)
        hand_x_temp = fltarr(npeak)
        hand_y_temp = fltarr(npeak)
        hand_fwhm_temp = fltarr(npeak)
        hand_sub_temp = lonarr(npeak)
        hand_minx_temp = fltarr(npeak)
        hand_maxx_temp = fltarr(npeak)
    ENDIF ELSE BEGIN
        handaps = 0
        hand_x_temp = 0
        hand_y_temp = 0
        hand_fwhm_temp = 0
        hand_sub_temp = 0
        hand_minx_temp = 0
        hand_maxx_temp = 0
    ENDELSE
    IF KEYWORD_SET(HAND_X) THEN BEGIN
       ;; Are there specified apertures for this slit?
       indslit = WHERE(handslits EQ slitid, nhand)
       
       IF nhand GT 0 THEN BEGIN
          ;; Add these aps to the xcen array
          xcen_hand = diag_matrix(bilinear(ximg, hand_x[indslit] $
                                           , hand_y[indslit])*nsamp)
          IF npeak GT 0 THEN xcen = [xcen, xcen_hand] $
          ELSE xcen = xcen_hand
          IF npeak GT 0 THEN ypeak = $
             [ypeak, interpolate(fluxconv, xcen_hand)] $
          ELSE ypeak = interpolate(fluxconv, xcen_hand)
          IF npeak GT 0 THEN handaps = [handaps, lonarr(nhand)+1] $
          ELSE handaps = lonarr(nhand)+1
          IF npeak GT 0 THEN hand_x_temp = [hand_x_temp, hand_x[indslit]] $
          ELSE hand_x_temp =  hand_x[indslit]
          IF npeak GT 0 THEN hand_y_temp = [hand_y_temp, hand_y[indslit]] $
          ELSE hand_y_temp = hand_y[indslit]
          IF npeak GT 0 THEN hand_fwhm_temp =  $
             [hand_fwhm_temp, hand_fwhm[indslit]] $
          ELSE  hand_fwhm_temp =  hand_fwhm[indslit]
          IF npeak GT 0 THEN hand_sub_temp =  $
             [hand_sub_temp, hand_sub[indslit]] $
          ELSE  hand_fwhm_temp =  hand_fwhm[indslit]
          IF npeak GT 0 THEN hand_minx_temp =  $
             [hand_minx_temp, hand_minx[indslit]] $
          ELSE  hand_minx_temp =  hand_minx[indslit]
          IF npeak GT 0 THEN hand_maxx_temp =  $
             [hand_maxx_temp, hand_maxx[indslit]] $
          ELSE  hand_maxx_temp =  hand_maxx[indslit]
          npeak = n_elements(xcen)
       ENDIF
    ENDIF

    ;stop
    ;;
    ;;  Go to next slit if no peaks are significant
    ;; 
    IF npeak EQ 0 THEN CONTINUE
    objstruct1 = long_obj_create(npeak, slitid = slitid, ny = ny)
    ;; Flag aps which were added by hand and set there x,y positions
    objstruct1.HAND_AP = handaps
    objstruct1.HAND_X = hand_x_temp
    objstruct1.HAND_Y = hand_y_temp
    objstruct1.HAND_FWHM = hand_fwhm_temp
    objstruct1.HAND_SUB = hand_sub_temp
    objstruct1.HAND_MINX = hand_minx_temp
    objstruct1.HAND_MAXX = hand_maxx_temp
    ;; Make a trace of the object, extrapolated to be in the same place
    ;; within the slit at all Y positions.

    FOR ipeak = 0L, npeak-1L DO BEGIN
        IF KEYWORD_SET(STDTRACE) THEN BEGIN
           ;; position of trace at middle of chip
           splog, 'Using input STDTRACE as crutch for object tracing'
            x_trace = interpol(STDTRACE, dindgen(n_elements(stdtrace)), ymid)
            shift = xx1[*, slitid-1] + (xcen[ipeak]/nsamp)*xsize[ymid] - x_trace
            objstruct1[ipeak].xpos = STDTRACE + shift
        ENDIF ELSE BEGIN
            ;; Below we set the first guess for the trace to be the
            ;; translated slit boundary, using the slit with the the
            ;; *largest* rms. The reason for doing this is that for
            ;; the LRIS blue chip, the slits are split at the chip
            ;; gap and set to be veritcal lines. To trace with the
            ;; proper curvature, the rms selection below will always
            ;; choose the curved trace of the other end of the slit.
            xx1_bar = total(xx1[*, slitid-1])/double(ny)
            xx2_bar = total(xx2[*, slitid-1])/double(ny)
            rms1 = sqrt(total((xx1[*, slitid-1]-xx1_bar)^2)/double(ny))
            rms2 = sqrt(total((xx2[*, slitid-1]-xx2_bar)^2)/double(ny))
            IF rms1 GT rms2 THEN objstruct1[ipeak].xpos = xx1[*, slitid-1] $
              + (xcen[ipeak] / nsamp) * xsize $
            ELSE  objstruct1[ipeak].xpos = xx2[*, slitid-1] $
              - (1.0 - xcen[ipeak]/nsamp)*xsize
         ENDELSE
        ;stop
        ;; Estimate the FWHM of this peak by literally looking for the
        ;; full width at the half-max (using linear interpolation)
        objstruct1[ipeak].xfracpos = xcen[ipeak]/nsamp
        objstruct1[ipeak].peakflux = ypeak[ipeak]
        yhalf = 0.5 * ypeak[ipeak]
        x0 = round(xcen[ipeak])
        if (x0 LT nsamp-1) then begin
            i2 = (where(fluxconv[x0:nsamp-1] LT yhalf))[0]
            IF i2 EQ -1 THEN xright = 0 $
            ELSE xright = $
              interpol(x0+[i2-1, i2], fluxconv[x0+i2-1:x0+i2], yhalf)
        endif else xright = 0
        if (x0 GT 0) then begin
           i1 = (reverse(where(fluxconv[0:x0] LT yhalf)))[0]
           IF i1 EQ -1 THEN xleft = 0 $
           ELSE begin 
              if i1 EQ (nsamp-1) then xleft = nsamp-1 else $
                 xleft = interpol([i1, i1+1], fluxconv[i1:i1+1], yhalf)
           endelse
        endif else xleft = 0
        fwhmmeas = 0.
         if (xleft NE 0 AND xright NE 0) then fwhmmeas = (xright - xleft) $
         else if (xleft NE 0) then fwhmmeas = 2*(xcen[ipeak] - xleft) $
         else if (xright NE 0) then fwhmmeas = 2*(xright - xcen[ipeak]) 

         if fwhmmeas NE 0 then objstruct1[ipeak].fwhm = $
           sqrt((fwhmmeas^2 - fwhm^2) > 4.) $
         else objstruct1[ipeak].fwhm = fwhm
      endfor
    ;stop
      ; Resort the traces by their X position
      objstruct1 = objstruct1[sort(objstruct1.xfracpos)]
      objstruct = keyword_set(objstruct) ? [objstruct, objstruct1] : objstruct1
      xtmp = (findgen(nsamp)+0.5) / nsamp
      qobj = xtmp*0.

      ii = where(objstruct.slitid EQ slitid, nobj)
      this_slit = where(slitmask EQ slitid, n_slit)
      FOR iobj = 0, nobj-1L DO BEGIN
         if skythresh gt 0 then objstruct[ii[iobj]].maskwidth = 3.*objstruct[ii[iobj]].fwhm * $
            (1. + 0.5*alog10((objstruct[ii[iobj]].peakflux/skythresh) > 1)) $
          else objstruct[ii[iobj]].maskwidth = 3.*objstruct[ii[iobj]].fwhm
          sep = abs(xtmp-objstruct[ii[iobj]].xfracpos)
          sep_inc = objstruct[ii[iobj]].maskwidth/nsamp 
          close = where(sep LE sep_inc, nc)
          if nc GT 0 then qobj[close] = qobj[close] + $
            objstruct[ii[iobj]].peakflux * $
            exp(-2.77*(sep[close]*nsamp)^2/(objstruct[ii[iobj]].fwhm)^2 > (-9.))
      ENDFOR

      if arg_present(objmask) then $
        objmask[this_slit] = objmask[this_slit] +  $
                   (interpol(qobj, xtmp, ximg[this_slit]) GT OBJTHRESH)

      ;----------
      ; Optionally find the best sky pixels
      if arg_present(skymask_new) then begin
          f = (qobj + fluxconv)[1:nsamp-2]
          sky_sort = sort(f)
          
          sky_ll = 2L
          xtmp2 = (findgen(nsamp+2)-0.5) / nsamp

;
;  use at least 1/4 pixels for sky 
;
          IF NOT KEYWORD_SET(SIGMA_SKY) THEN SIGMA_SKY = 7.0D
          sky_uu = min(where(f[sky_sort] GE sigma_sky*skythresh)) > $
                   (long(nsamp * 0.25)+2L)  

         oksky = xtmp2*0.
         oksky[sky_sort[sky_ll:sky_uu]+2] = 1.

         skymask[this_slit] = skymask[this_slit] * $
                   (interpol(oksky, xtmp2, ximg[this_slit]) GT 0.5)
     endif
 endfor


   if (arg_present(skymask)) then begin
      nskypix = long(total(skymask))
      slit_pix = WHERE(slitmask NE 0, n_slitpix)
      splog, 'Fraction of image for sky = ', float(nskypix)/double(n_slitpix)
   endif

   if NOT keyword_set(objstruct) then begin
      splog, 'Number of objects =  0'
      return, 0
   endif

   ;;
   ;; Tweak traces and flag bad fits
   ;;
   IF KEYWORD_SET(VERBOSE) THEN splog, 'Tweaking the object positions'
   
   ;; First due a trace_crude to lock onto the objects. Fit this
   ;; crude trace and then use that as the input trace for flux and 
   ;; then gaussian weighted centering. 
                                ;fwhm2 = djs_median(objstruct.FWHM)
   IF KEYWORD_SET(CRUDE) THEN BEGIN
      xpos_slit = objstruct.xpos
      ypos = objstruct.ypos
      xstart = djs_median(xpos_slit[(ymid-10L):(ymid+10L), *], 1)
      xcrude = trace_crude(image*mask, xstart = xstart, ystart = ymid $
                           , radius = fwhm, nave = nave $
                           , maxshifte = maxshifte $
                           , maxshift0 = maxshift0, xerr = xerr_crude)
      xy2traceset, ypos, xcrude, ss_guess $
                   , invvar = 1./(xerr_crude^2 + (xerr_crude EQ 0.0)) $
                   , func = func_crude, ncoeff = ncoeff_crude, yfit = xpos0 $
                   , silent = silent 
      xdiff_med = djs_median(xpos_slit-xpos0, 1)
;  If the median difference between slit trace and crude trace is too big
;  then just use the slit trace
      FOR kobj = 0L, nobj-1L DO $
         IF abs(xdiff_med[kobj]) GT CRUDE_TOL $
         THEN xpos0[*, kobj] = xpos_slit[*, kobj]
   ENDIF ELSE BEGIN
      xpos0 = objstruct.xpos
      ypos  = objstruct.ypos
   ENDELSE
;   stop
   ;; Fixed bug where invvar's were being used. a single invvar=0 pixel 
   ;; returns as the fweigthed (or gweighted) trace to initial guess. 
   niter = 12L
   xfit1 = xpos0
   FOR i = 1L, niter DO BEGIN
      IF i LT niter/3 THEN fwhm_now = 1.3*fwhm $
      ELSE IF (i GE niter/3) AND (i LT 2*niter/3) THEN fwhm_now = 1.1*fwhm $
      ELSE fwhm_now = fwhm
      xpos1 = trace_fweight(image*mask, xfit1, ypos, radius = fwhm_now)
      ;; set these threhsolds to avoid poorly behaved fits
      xpos1 = xpos1 > (-0.2*nx)
      xpos1 = xpos1 <  1.2*nx
      ;; Use invvar to deal with truncated slits (e.g. ESI)  :: 
      ;; JXP 19 Feb 2015
      if keyword_set(CLIP_SLIT) then $
         xinvvar = long_obj_clipslit( thisimg, clip_slit, xpos1, ypos, objstruct)
      ;; Fit
      xy2traceset, ypos, xpos1, pos_set1, ncoeff = ncoeff, yfit = xfit1 $
                   , maxdev = 5.0, /silent, invvar=xinvvar
   ENDFOR

   xfit2 = xfit1
   FOR i = 1L, niter DO BEGIN
      xpos2 = trace_gweight(image*mask, xfit2, ypos, sigma = fwhm/2.3548)
      ;; set these threhsolds to avoid poorly behaved fits
      xpos2 = xpos2 > (-0.2*nx)
      xpos2 = xpos2 <  1.2*nx
      ;; Use invvar to deal with truncated slits (e.g. ESI)  :: 
      ;; JXP 19 Feb 2015
      if keyword_set(CLIP_SLIT) then $
         xinvvar = long_obj_clipslit( thisimg, clip_slit, xpos2, ypos, objstruct)
      ;;
      xy2traceset, ypos, xpos2, pos_set2, ncoeff = ncoeff, yfit = xfit2 $
                   , maxdev = 5.0, /silent, invvar=xinvvar
   ENDFOR
;   xposivar = 0.1/(xerr2^2 + 0.02^2) * (xerr2 GT 0 AND xerr2 LT 2.0)
   ;; More robust to not include the formal centroiding errors
   ;; below...
   ;; Use invvar to deal with truncated slits (e.g. ESI)  :: 
   ;; JXP 19 Feb 2015
   if keyword_set(CLIP_SLIT) then $
      xinvvar = long_obj_clipslit( thisimg, clip_slit, xpos2, ypos, objstruct)
   xy2traceset, ypos, xpos2, pos_set, ncoeff = ncoeff, yfit = xfit $
                , maxdev = 5.0, /silent, invvar=xinvvar
   objstruct.xpos = xfit

   ;; Parse out bad objects -- JXP March 2009
   nobj = n_elements(objstruct)
   msk = bytarr(nobj)
   for kk=0L,nobj-1 do begin
       mx = max(objstruct[kk].xpos)
       if mx GT 1e-4 then msk[kk] = 1B
   endfor
   objstruct = objstruct[where(msk,nobj)]
   objstruct.objid = lindgen(nobj) + 1L

   ;; Now for the hand_aps go through and assign the trace to be 
   ;; the brightest object on the slit, or use the slit boundary if there
   ;; is no object

   IF KEYWORD_SET(HAND_X) THEN BEGIN
      FOR slitid=1L, max(slitmask) do begin
           handinds = WHERE(objstruct.HAND_AP EQ 1 $
                           AND objstruct.SLITID EQ slitid, nhand)
           reguinds = WHERE(objstruct.HAND_AP EQ 0 $
                           AND objstruct.SLITID EQ slitid, nreg)
           IF nhand GT 0 THEN BEGIN
               IF nreg GT 0 THEN BEGIN ; find the brightest
                   brightest = max(objstruct[reguinds].PEAKFLUX, imax)
                   x_0 = interpol(objstruct[reguinds[imax]].XPOS $
                                  , objstruct[reguinds[imax]].YPOS $
                                  , objstruct[handinds].HAND_Y)
                   shift = objstruct[handinds].HAND_X - x_0
                   objstruct[handinds].xpos =  $
                     objstruct[reguinds[imax]].XPOS # replicate(1.0, nhand) + $
                     replicate(1.0, ny) # shift
               ENDIF ELSE BEGIN 
                  ;; instead use either the slit boundary or input trace
                  xsize = xx2[*, slitid-1] - xx1[*, slitid-1]
                  IF KEYWORD_SET(STDTRACE) THEN BEGIN
                     ;; position of trace at middle of chip
                     x_0 = interpol(STDTRACE $
                                    , dindgen(n_elements(stdtrace)) $
                                    , objstruct[handinds].HAND_Y)
                     shift = objstruct[handinds].HAND_X - x_0
                     objstruct[handinds].xpos = $
                        STDTRACE # replicate(1.0, nhand) $
                        + replicate(1.0, ny) # shift
                  ENDIF ELSE BEGIN
                     xx1_bar = total(xx1[*, slitid-1])/double(ny)
                     xx2_bar = total(xx2[*, slitid-1])/double(ny)
                     rms1 = sqrt(total((xx1[*, slitid-1]-xx1_bar)^2)/ $
                                 double(ny))
                     rms2 = sqrt(total((xx2[*, slitid-1]-xx2_bar)^2)/ $
                                 double(ny))
                     IF rms1 GT rms2 THEN objstruct[handinds].xpos = $
                        xx1[*, slitid-1] # replicate(1.0, nhand) +   $
                        xsize # objstruct[handinds].XFRACPOS $
                     ELSE objstruct[handinds].xpos = $
                        xx2[*, slitid-1] # replicate(1.0, nhand) -   $
                        xsize # (1.0 - objstruct[handinds].XFRACPOS)
                  ENDELSE
               ENDELSE
               ;; loop over the regular aps and exclude any which are
               ;; within a FWHM of the hand_aps
               IF nreg GT 0 THEN BEGIN
                  msk = lonarr(n_elements(objstruct)) + 1L
                  if nreg EQ 1 then reg_fwhm = objstruct[reguinds].FWHM else $
                     reg_fwhm = median(objstruct[reguinds].FWHM)
                  reg_xpos = djs_median(objstruct[reguinds].xpos, 1)
                  han_xpos = djs_median(objstruct[handinds].xpos, 1)
                  FOR ii = 0L, nreg-1L DO BEGIN
                     iclose = WHERE(abs(reg_xpos[ii] - han_xpos) LE 0.6*reg_fwhm, nclose)
                     IF nclose GT 0 THEN msk[reguinds[ii]] = 0L
                  ENDFOR
                  objstruct = objstruct[where(msk, nobj)]
                  objstruct.objid = lindgen(nobj) + 1L
               ENDIF
            ENDIF
        ENDFOR
      FOR kk = 0L, n_elements(objstruct)-1L DO BEGIN
         IF objstruct[kk].HAND_AP EQ 1 THEN BEGIN
            IF objstruct[kk].HAND_MINX GT 0 OR $
               objstruct[kk].HAND_MAXX GT 0 THEN BEGIN
               max_sep = $
                  max([objstruct[kk].HAND_MINX, objstruct[kk].HAND_MAXX])
               objstruct[kk].MASKWIDTH = 1.5*max_sep
               objstruct[kk].FWHM = objstruct[kk].HAND_FWHM
            ENDIF ELSE IF objstruct[kk].HAND_FWHM NE 0.0 THEN BEGIN
               objstruct[kk].FWHM = objstruct[kk].HAND_FWHM
               objstruct[kk].MASKWIDTH = 3.0d*objstruct[kk].HAND_FWHM
            ENDIF
         ENDIF
      ENDFOR
   ENDIF

   median_fwhm = djs_median(objstruct.FWHM)
   skymask_new = slitmask GT 0
   nobj = n_elements(objstruct)
   left_ind  = objstruct.xpos - $
     replicate(median_fwhm, nobj) ## replicate(1.0D, ny)/2.0D
   right_ind = objstruct.xpos + $
     replicate(median_fwhm, nobj) ## replicate(1.0D, ny)/2.0D
   FOR iobj = 0L, n_elements(objstruct) - 1L DO BEGIN
       FOR j = 0L, ny-1L DO BEGIN 
           xmin = (floor(left_ind[j, iobj]) >  0) <  (nx-2L)
           xmax = (ceil(right_ind[j, iobj]) <  (nx-1L)) > 0 
           skymask_new[xmin:xmax, j] = 0L
       ENDFOR
    ENDFOR

   nobj = n_elements(objstruct)
   objstruct.objid = lindgen(nobj) + 1L

   IF NOT KEYWORD_SET(SILENT) THEN BEGIN
      splog, 'Number of objects = ', n_elements(objstruct)
      for i = 0L, n_elements(objstruct)-1L do $
         splog, 'objid ', i+1, ', slitid ', objstruct[i].slitid 
      ;;, ', tweak: ', $
      ;;          median(abs(xfit[*, i]-objstruct[i].xpos)), ' pixels'
   ENDIF

   return, objstruct
end

;------------------------------------------------------------------------------
