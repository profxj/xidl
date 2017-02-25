;+
; NAME:
;  long_combspec 
;  Version 1.1
;
; PURPOSE:
;  Combine the spectra from multiple exposures taken through the same
;  mask or longslit.
;
; CALLING SEQUENCE:
;  LONG_COMBSPEC,  influx, inivar, loglam, [ newloglam = , newflux = $
;                  , newivar = , newmask = , iref = , SIGREJ = ]
;
; INPUTS:
;   influx     -- Flux array with dimension [nspec,nimgs]
;   inivar     -- Inverse variance array with dimension [nspec,nimgs]
;   loglam     -- Log10 wavelength array. Can be either dimension 
;                 [nspec] or dimension [nspec,nimgs]. If the loglam
;                 array is passed as a 2d array then the code will perform
;                 spline interpolation onto the common wavelength grid
;                 of index iref in this array. 
;
; OPTIONAL INPUTS:
;   iref       -- Index for reference of scaled spectra (and wavelength grid)
;                 if 2-d wavelength array is passed. 
;   SIGREJ     -- Parameter for rejection and masking. Default is SIGREJ=3
;
;   /CHECK     -- Plot the scaled data against the reference. 
;   newloglam  -- One can also input a 'final' wavelength array that
;                 (presumably) differs from the input arrays
;   INPUT_WEIGHTS -- User specified weights for the co-addition. This
;                    has to be a vector of floats of size [nimgs]
; OUTPUTS:
; 
; OPTIONAL OUTPUTS:
;   newloglam  -- New wavelength array. Same as loglam if loglam is a 1d 
;                 array. Otherwise newloglam=loglam[*,iref]
;   newflux    -- Combined flux
;   newivar    -- Combined inverse variance
;   newmask    -- Mask for combined spectra. 
;   arrflux    -- The [nspec,nimgs] array holding scaled flux vectors
;   arrivar    -- The [nspec,nimgs] array holding scaled ivar vectors
;
; COMMENTS:
;   This code will interpolate spectra with different wavelength vectors
;   onto a common grid, or  a single wavelength vector can be passed. It 
;   uses polynomial scaling to the first reference spectrum, performs
;   sigma clipping to mask outlier pixels, and then performs a SNR waited
;   average of the spectra.
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;  
;
; REVISION HISTORY:
;   18-May-2007  Created by JFH
;-
;------------------------------------------------------------------------------
;    exten=  --  Extension for the spectral file (default = 4)

PRO RESID_GAUSS_PLOT, chi, one_sigma

max = 6.0
min = -6.0
n_bins = 50
n_tot = n_elements(chi)
binsize = (max-min)/float(n_bins)
bins_histo = min + findgen(n_bins)*binsize + binsize/2.0 
histo = float(HISTOGRAM(chi, BINSIZE = binsize $
                        , MIN = min, NBINS = n_bins))
n_histo = total(histo)
norm = float(n_histo)*binsize
histo = histo/norm
xvals = -10.0 + 0.02*findgen(1001)
ygauss = gauss1(xvals, [0.0, 1.0, 1.0])
ygauss_new = gauss1(xvals, [0.0, one_sigma, 1.0])

x_splot, bins_histo, histo, psym1 = 10, psym2 = 3 $
         , xtwo = xvals, ytwo = ygauss, xthr = xvals, ythr = ygauss_new $
         , psym3 = 3, title = 'Residual Distribution. New sigma=' $
         + strcompress(string(one_sigma, format = '(F7.2)'), /rem) $
         , XMNX = [-6.0, 6.0], YMNX = [0.0D, 0.7D], /block
RETURN
END

PRO LONG_COMBSPEC, influx, inivar, loglam $
                   , INSKY = INSKY, INNIVAR = INNIVAR $
                   , newloglam = newloglam, newflux = newflux $
                   , newivar = newivar, newmask = newmask $
                   , newnivar = newnivar, newsky = newsky $
                   , iref = iref1, SIGREJ = SIGREJ1, CHECK = CHECK $
                   , NOREJ = NOREJ, NOSHIFT = NOSHIFT, SHIFT_SN = SHIFT_SN1 $
                   , NOSHARP = NOSHARP $
                   , MASKLAM = MASKLAM1, LOGLAMMASK = LOGLAMMASK $
                   , LAM_MASK_MIN = LAM_MASK_MIN1, LAM_MASK_MAX = LAM_MASK_MAX1 $
                   , MEDSCALE = MEDSCALE, SN_MIN_MEDSCALE = SN_MIN_MEDSCALE1 $
                   , SN_MAX_MEDSCALE = SN_MAX_MEDSCALE1 $
                   , DEBUG = DEBUG $
                   , SN2 = SN2_WEIGHTS, YMULT = YMULT $
                   , HAND_SCALE = HAND_SCALE, EQUAL_WEIGHT = EQUAL_WEIGHT $
                   , DV = dv, MEAN_SN2 = MEAN_SN2, IN_NPOLY = IN_NPOLY $
                   , FMIN = fmin, FMAX = fmax, MXSHIFT = MXSHIFT1, XCORR_SKY = XCORR_SKY $
                   , POLY_RATIO_SN = POLY_RATIO_SN $
                   , NO_OFFSET = NO_OFFSET, USE_AVG_SN_WEIGHTS = USE_AVG_SN_WEIGHTS $
                   , INPUT_WEIGHTS = INPUT_WEIGHTS

IF NOT KEYWORD_SET(SIGREJ1) THEN sigrej = 3.0D $
ELSE sigrej = sigrej1
IF KEYWORD_SET(NBACK1) THEN nback = nback1 ELSE nback = 0
IF NOT KEYWORD_SET(INSKY) THEN INSKY = 0.0*influx
IF NOT KEYWORD_SET(INNIVAR) THEN INNIVAR = 0.0*influx
IF N_ELEMENTS(SHIFT_SN1) GT 0 THEN SHIFT_SN = SHIFT_SN1 ELSE SHIFT_SN = 7.0d
IF KEYWORD_SET(SN_MIN_MEDSCALE1) THEN SN_MIN_MEDSCALE = SN_MIN_MEDSCALE1 ELSE SN_MIN_MEDSCALE = 0.5
IF KEYWORD_SET(SN_MAX_MEDSCALE1) THEN SN_MAX_MEDSCALE = SN_MAX_MEDSCALE1 ELSE SN_MAX_MEDSCALE = 2.0
IF NOT KEYWORD_SET(SCL_MIN) THEN SCL_MIN = 0.1
IF NOT KEYWORD_SET(SCL_MAX) THEN SCL_MAX = 100.0


if not keyword_set(DV) then dv = 10000.0 ;; med width in km/s
bkspace = (dv/3.0d5)/alog(10.0d)

sigrej_int = 3.0                
; internal sigrej for djs_iterstat and djs_avsigclip calls
dims = size(influx, /dim)
IF n_elements(dims) EQ 1 THEN BEGIN
    nimgs = 1
    nsz = dims[0]
    iref = 0
ENDIF ELSE BEGIN
    nimgs = dims[1]
    nsz = dims[0]
    IF n_elements(iref1) GT 0 THEN iref = iref1 $
    ELSE iref = 0
ENDELSE

sizelam = size(loglam)
dims_lam = sizelam[0]

IF nimgs EQ 1 THEN BEGIN
    newloglam = loglam[*, iref]
    newflux = influx[*, iref]
    newivar = inivar[*, iref]
    newmask = (inivar[*, iref] GT 0)
    newsky  =  insky[*, iref]
    newnivar = innivar[*, iref]
    RETURN
 ENDIF

IF nimgs GT 1 AND dims_lam EQ 1 THEN $
   loglam = loglam # replicate(1.0d, nimgs)

IF KEYWORD_SET(HAND_SCALE) AND n_elements(HAND_SCALE) NE nimgs THEN $
   message $
   , 'ERROR: For hand scaling number of imgs must equal n_elements(HAND_SCALE)'


IF NOT KEYWORD_SET(NEWLOGLAM) THEN BEGIN
;; construct new wavelength grid by concatenating on new wavelengths
   IF dims_lam GT 1 THEN BEGIN
      newloglam = loglam[*, iref]
      FOR j = 0L, nimgs-1 DO BEGIN
         IF j EQ iref THEN CONTINUE
         nspec = n_elements(newloglam)
         dloglam_0 = (newloglam[1]-newloglam[0])
         dloglam_n =  (newloglam[nspec-1L] - newloglam[nspec-2L])
         IF (newloglam[0] - loglam[0, j]) GT dloglam_0 THEN BEGIN
            min1 = min(abs(loglam[*, j] - newloglam[0] - dloglam_0), kmin)
            newloglam = [loglam[0:kmin, j], newloglam]
         ENDIF
         nspec = n_elements(newloglam)
         IF (loglam[nsz-1L, j] - newloglam[nspec-1L]) GT dloglam_n THEN BEGIN
            min1 = min(abs(loglam[*, j] - newloglam[nspec-1L] - dloglam_n) $
                       , kmin)
            newloglam = [newloglam, loglam[kmin:*, j]]
         ENDIF
      ENDFOR
   ENDIF ELSE  newloglam = loglam[*, iref]
ENDIF


nspec = n_elements(newloglam)
IF KEYWORD_SET(MASKLAM1) AND KEYWORD_SET(LOGLAMMASK) THEN BEGIN
    masklam2 = interpol(double(masklam1), loglammask, newloglam)
    masklam = masklam2 GT 0.9D
ENDIF ELSE IF KEYWORD_SET(MASKLAM1) OR KEYWORD_SET(LOGLAMMASK) THEN $
  message, 'Must specify both masklam and loglammask together' $
ELSE BEGIN
   masklam = lonarr(nspec)
   IF NOT KEYWORD_SET(LAM_MASK_MIN1) THEN LAM_MASK_MIN = 0.0 ELSE LAM_MASK_MIN = LAM_MASK_MIN1 
   IF NOT KEYWORD_SET(LAM_MASK_MAX1) THEN LAM_MASK_MAX = 1.0d10 ELSE LAM_MASK_MAX = LAM_MASK_MAX1 
   igd = WHERE(10.0d^newloglam GT LAM_MASK_MIN AND 10.0d^newloglam LT LAM_MASK_MAX, ngd)
   IF ngd EQ 0 THEN message, 'No wavelengths in specified range LAM_MASK_MIN and LAM_MASK_MAX' $
   ELSE masklam[igd] = 1
ENDELSE

arrflux   = dblarr(nspec, nimgs)
arrivar   = dblarr(nspec, nimgs)
arrmask   = dblarr(nspec, nimgs)
outmask   = dblarr(nspec, nimgs)
arrsky    = dblarr(nspec, nimgs)
arrnivar  = dblarr(nspec, nimgs)
sclflux   = dblarr(nspec, nimgs)
sclivar   = dblarr(nspec, nimgs)
sclmask   = intarr(nspec, nimgs) ;; JXP 03 Oct 2014
sclsky    = dblarr(nspec, nimgs)
sclnivar  = dblarr(nspec, nimgs)
ymult     = dblarr(nspec, nimgs)
sn2_weights = dblarr(nspec, nimgs)
sn_weights = dblarr(nspec, nimgs)
mean_sn2 = dblarr(nimgs)

; map the loglam to interval [0,1] for polynomial scaling
diff_loglam = max(newloglam)-min(newloglam)
xvector = (newloglam-min(newloglam))/diff_loglam
;; Rebin spectra onto the same wavelength grid
FOR j = 0L, nimgs-1L DO BEGIN
    ;; Do we need to rebin the wavelengths??
    IF dims_lam EQ 1 THEN BEGIN
        sclflux[*, j] = influx[*, j]
        sclivar[*, j] = inivar[*, j]
        mask1 = double((inivar[*, j] GT 0.0) $
                       AND (abs(influx[*, j]) LE 1.0d6) $
                       AND (finite(influx[*, j]) EQ 1)  $
                       AND (finite(inivar[*, j]) EQ 1)  $
                       AND (inivar[*, j] LT 1.0d8))
        IF NOT KEYWORD_SET(NOSHARP) THEN BEGIN
            spix = WHERE(mask1)
            sharpchi =  ((influx[*, j]-djs_median(influx[*, j] $
                                                     , width = 3 $
                                                     , bound = 'reflect'))  $
                         *sqrt(inivar[*, j] >  0.0))
            djs_iterstat, sharpchi[spix], mean = mean, median = median $
                          , sigma = sigma
            ;; filter out small values in case the mode is dominating the stats
            spix = WHERE(mask1 AND abs(sharpchi) GT ((0.2*sigma) > 0.1D))
            djs_iterstat, sharpchi[spix], mean = mean, median = median $
                          , sigma = sigma
            sharpmask = sharpchi LT ((median >  mean) + 10.0*sigma)
            sharpmask = long_grow_mask(sharpmask, 1)
            mask = mask1 AND sharpmask ; 0=bad, 1=good
        ENDIF ELSE mask = mask1
        sclmask[*, j] = (mask EQ 0) ; 0=good, 1=bad
        sclsky[*, j]  = insky[*, j]
        sclnivar[*, j] = innivar[*, j]
    ENDIF ELSE BEGIN
;;;;;;;Change 2: exclude the points we have already set to zero
        in_lam = WHERE(newloglam GE min(loglam[where(loglam[*,j] gt 0.), j]) AND $
                       newloglam LE max(loglam[*, j]) $
                       , COMPLEMENT = out_lam, NCOMPLEMENT = nout)
        ;; interpolate the mask onto the new wavelength grid
        ;; and mask any pixels not covered by this spectrum
        mask1 = double((inivar[*, j] GT 0.0) $
                       AND (abs(influx[*, j]) LE 1.0d6) $
                       AND (finite(influx[*, j]) EQ 1)  $
                       AND (finite(inivar[*, j]) EQ 1)  $
                       AND (inivar[*, j] LT 1.0d8))
        IF NOT KEYWORD_SET(NOSHARP) THEN BEGIN
            spix = WHERE(mask1)
            sharpchi =  ((influx[*, j]-djs_median(influx[*, j], width = 3 $
                                                  , bound = 'reflect'))  $
                         *sqrt(inivar[*, j] >  0.0))            
            djs_iterstat, sharpchi[spix], mean = mean, median = median $
                          , sigma = sigma
            spix = WHERE(mask1 AND abs(sharpchi) GT  ((0.2*sigma) > 0.1D))
            djs_iterstat, sharpchi[spix], mean = mean, median = median $
                          , sigma = sigma
            sharpmask = sharpchi LT (median + 10.0*sigma)
            sharpmask = long_grow_mask(sharpmask, 1) ;; grow mask by 1
            mask = mask1 AND sharpmask
        ENDIF ELSE mask = mask1
        sclmask1 = interpol(double(mask), loglam[*, j], newloglam[in_lam])
        sclmask[in_lam, j]  = (sclmask1 LE 0.5D)
        IF nout GT 0 THEN sclmask[out_lam, j] = 1
        ;; interpolate over the masked pixels for wave rebinning,
        ;; then make sure bad pixels are masked. 
        ;; was using quadratic interpolation before but this is unstable 
        ;; near sky lines so I changed to linear interpolation which 
        ;; looks more stable than all other alternatives. 

        ;; Deal with zero values (JXP 05 Oct 2014)
        gdp = where(loglam[*,j] GT 0.)
        ivar1    = djs_maskinterp(inivar[*, j], mask EQ 0)

        sclivar[in_lam, j]  = interpol(ivar1[gdp], loglam[gdp, j] $
                                       , newloglam[in_lam])* $
          (sclmask[in_lam, j] EQ 0) 
        ;; this absolute value line approximately deals with cases 
        ;; where the interpolated ivar crosses zero because of interpolation
        ;; error. 
        sclivar[in_lam, j] = abs(sclivar[in_lam, j])
        sclmask[in_lam, j] = sclmask[in_lam, j] OR (sclivar[in_lam, j] EQ 0.0)
        flux1    = djs_maskinterp(influx[*, j], mask EQ 0)
        sclflux[in_lam, j]  = interpol(flux1[gdp], loglam[gdp, j], newloglam[in_lam]) $
                              *(sclmask[in_lam, j] EQ 0)
        sky1     = djs_maskinterp(insky[*, j], mask EQ 0)
        sclsky[in_lam, j]   = interpol(sky1[gdp], loglam[gdp, j], newloglam[in_lam]) $
          *(sclmask[in_lam, j] EQ 0)
        nivar1   = djs_maskinterp(innivar[*, j], mask EQ 0)
        sclnivar[in_lam, j] = interpol(nivar1[gdp], loglam[gdp, j], $
                                       newloglam[in_lam])*(sclmask[in_lam, j] EQ 0)
    ENDELSE
 ENDFOR

FOR j = 0L, nimgs-1L DO BEGIN
   ;; For the mean S/N computation only consider wavelengths in the masklam
   gdp = where(loglam[*,j] GT 0.)
   in_lam_mean = WHERE(newloglam GE min(loglam[gdp, j]) AND $
                       newloglam LE max(loglam[*, j]) AND $
                       masklam, nin $
                       , COMPLEMENT = out_lam, NCOMPLEMENT = nout)
   sn2_mean = (sclflux[in_lam_mean, j]^2*sclivar[in_lam_mean, j] >  0.0)
   ind_nonzero = where(sclivar[in_lam_mean, j] GT 0 AND sn2_mean GT 0, nnozero)
   IF nnozero GT 0 THEN djs_iterstat, sn2_mean[ind_nonzero], mean = mean_sn21 $
                                      , sigrej = sigrej_int, sigma = sigma_sn $
   ELSE mean_sn21 = 0.0
   mean_sn2[j] = mean_sn21
   ;; Now for the weight functions use the entire wavelength array
   in_lam = WHERE(newloglam GE min(loglam[*, j]) AND $
                  newloglam LE max(loglam[*, j]), nin $
                  , COMPLEMENT = out_lam, NCOMPLEMENT = nout)
   sn2 = (sclflux[in_lam, j]^2*sclivar[in_lam, j] >  0.0)
   bsp_mask =  (sclmask[in_lam, j] EQ 0) AND (sclivar[in_lam, j] GT 0) $
               AND (sn2 GT 0) AND (sclflux[in_lam, j] GT 0.0)
   ibsp_good = WHERE(bsp_mask, nbsp)
   IF nbsp EQ 0 THEN message, 'LONG_COMBSPEC: Found no good pixels in exposure #' + $ ; '
                              strcompress(string(j), /rem) + $
                              '. There is likely something wrong with this data'
   med_width = round(nbsp/((max(newloglam[ibsp_good]) $
                            - min(newloglam[ibsp_good]))/bkspace))
   ;;med_width = round(nbsp/((max(newloglam[in_lam]) $
   ;;                         - min(newloglam[in_lam]))/bkspace))
   sn2_med = fltarr(nin)

   sn2_med1 = djs_median(sn2[ibsp_good], width = med_width $
                         , boundary = 'reflect')
   sn2_med2 = interpol(sn2_med1, newloglam[in_lam[ibsp_good]] $
                       , newloglam[in_lam])
   
   sig_res = med_width/10L > 3L
   nhalf =  long(sig_res)*4L
   xkern = dindgen(2*nhalf+1)-nhalf
   kernel = gauss1(xkern, [0.0, sig_res, 1.0])
   sn2_med = convol(sn2_med2, kernel, /edge_truncate)
   sn2_weights[in_lam, j] = sn2_med >  1.0D ;; Do we really want to set a floor here???
   print, 'sqrt(S/N^2): ', j, sqrt(mean_sn2[j])


   ;; S/N calculation (not S/N^2!)
   sn = (sclflux[in_lam, j]*sqrt(sclivar[in_lam, j]))
   bsp_mask_sn =  (sclmask[in_lam, j] EQ 0) AND (sclivar[in_lam, j] GT 0)
   ibsp_good_sn = WHERE(bsp_mask_sn, nbsp_sn)
   IF nbsp_sn EQ 0 THEN message, 'LONG_COMBSPEC: Found no good pixels in exposure #' + $
                                 strcompress(string(j), /rem) + $ ;'
                                 '. There is likely something wrong with this data'
   med_width_sn = round(nbsp_sn/((max(newloglam[ibsp_good]) $
                                  - min(newloglam[ibsp_good]))/bkspace))
   sn_med = fltarr(nin)

   
   sn_med1 = djs_median(sn[ibsp_good_sn], width = med_width $
                        , boundary = 'reflect')
   sn_med2 = interpol(sn_med1, newloglam[in_lam[ibsp_good_sn]] $
                      , newloglam[in_lam])
   sig_res_sn = med_width_sn/10L > 3L
   nhalf_sn =  long(sig_res_sn)*4L
   xkern_sn = dindgen(2*nhalf_sn+1)-nhalf_sn
   kernel_sn = gauss1(xkern_sn, [0.0, sig_res_sn, 1.0])
   sn_med = convol(sn_med2, kernel_sn, /edge_truncate)
   sn_weights[in_lam, j] = sn_med  ;;>  1.0D


   ;; TESTING
   ;;junk = where(finite(sn2_med) EQ 0, njunk)
   ;;IF njunk GT 0 THEN stop
   ;stop
   IF KEYWORD_SET(DEBUG) THEN BEGIN
      x_splot, 10.0d^newloglam[in_lam], sn2 $
               ,  psym1 = 1, /block $ 
               , title = '(S/N)^2 fit for' + strcompress(string(j), /rem) $
               + 'th exposure', XMNX = 10.0D^[min(newloglam[in_lam]) $
                                              , max(newloglam[in_lam])] $
               , YMNX = [-1.0D, 1.2*max(sn2_weights[in_lam, j])] $
               , xtwo = 10.0d^newloglam, ytwo = sn2_weights[*, j] $
               , psym2 = 0
      ;;x_splot, 10.0d^newloglam[in_lam], sn $
      ;;         ,  psym1 = 1, /block $ 
      ;;         , title = '(S/N) fit for' + strcompress(string(j), /rem) $
      ;;         + 'th exposure', XMNX = 10.0D^[min(newloglam[in_lam]) $
      ;;                                        , max(newloglam[in_lam])] $
      ;;         , YMNX = [-1.0D, 1.2*max(sn_weights[in_lam, j])] $
      ;;         , xtwo = 10.0d^newloglam, ytwo = sn_weights[*, j] $
      ;;         , psym2 = 0
    ENDIF
ENDFOR

avg_sn = sqrt(total(mean_sn2)/double(nimgs))
splog, 'Average S/N for the image stack =' + string(avg_sn)
IF avg_sn LE 4.0D OR KEYWORD_SET(USE_AVG_SN_WEIGHTS) THEN BEGIN
   sn2_weights = (replicate(1.0D, nspec) # mean_sn2)
   splog, 'Low S/N regime. Using one (S/N)^2 weight per spectrum'
ENDIF

;; Weight all pixels and all exposures equally
IF KEYWORD_SET(EQUAL_WEIGHT) THEN sn2_weights[*] = 1.0d

;; User specifies the weights 
IF KEYWORD_SET(INPUT_WEIGHTS) THEN BEGIN
   sn2_weights = (replicate(1.0D, nspec) # mean_sn2)
   splog, 'Using user specified weights for this image stack'
ENDIF

;; For the reference flux, avsigclip the data to obtain a high SNR
;; average. Scale all spectra to the average counts. The reference
;; refivar is only used for the polynomial fitting
IF n_elements(IREF1) GT 0 THEN BEGIN
    refflux = sclflux[*, iref]
    refivar = sclivar[*, iref]
    refmask = refivar GT 0.0
    refsky =  sclsky[*, iref]
ENDIF ELSE BEGIN
    refflux = djs_avsigclip(sclflux, 2, inmask = sclmask $
                            , outmask = outmask1, sigrej = sigrej_int)
; bad everywhere?
    refmask = total(outmask1, 2) NE nimgs ; 0=bad,1=good
    sig2 = 1.0/(sclivar + (sclivar LE 0.0))
    nused = total(outmask1 EQ 0, 2)
    newsig2 = total(sig2*(outmask1 EQ 0), 2)/(nused^2 + (nused EQ 0))
    refivar = refmask/(newsig2 + (newsig2 LE 0.0))
    refsky = djs_avsigclip(sclsky, 2, inmask = sclmask $
                           , outmask = outmask2, sigrej = sigrej_int)
    
;; sharpness filter the combined spectrum
ENDELSE
goodpix = WHERE(refmask AND masklam, ngood)
IF ngood EQ 0 THEN goodpix = lindgen(nspec) ;; prevent crash
ref_smth = djs_median(refflux[goodpix], width = 5, bound = 'reflect')
ref_min = -0.1D                 ;min(ref_smth) > (-10.0)
ref_max = 1.5*max(ref_smth) < 1.0d6


IF NOT KEYWORD_SET(NOSHARP) THEN BEGIN
    sharpchi =  ((refflux-djs_median(refflux, width = 3 $
                                     , bound = 'reflect'))*sqrt(refivar >  0.0))
    spix = where(refmask)
    djs_iterstat, sharpchi[spix], mean = mean, median = median, sigma = sigma
    spix = where(refmask AND abs(sharpchi) GT ((0.2*sigma) > 0.1D))
    djs_iterstat, sharpchi[spix], mean = mean, median = median, sigma = sigma
    sharpmask = sharpchi LT (median + 10.0*sigma)
    sharpmask = long_grow_mask(sharpmask, 3)
    refmask = refmask AND sharpmask
ENDIF

refivar = refivar*refmask
goodpix = WHERE(refivar GT 0.0 AND finite(refflux) AND finite(refivar) $
                AND refmask EQ 1 AND refivar LT 1.0d8) 
djs_iterstat, refflux[goodpix], sigrej = sigrej_int, median = med_ref $
              , invvar = refivar[goodpix], mask = mask, sigma = sig_ref1
goodpix = WHERE(refivar GT 0.0 AND finite(refflux) AND finite(refivar) $
                AND refmask EQ 1 AND refflux GT 0.5*abs(med_ref) $
                AND refivar LT 1.0d8) 
djs_iterstat, refflux[goodpix], sigrej = sigrej_int, mean = med_ref $
              , invvar = refivar[goodpix], mask = mask, sigma = sig_ref1

IF KEYWORD_SET(MXSHIFT1) THEN MXSHIFT = MXSHIFT1 ELSE MXSHIFT = 2
nsamp = 50
step = lindgen(2*MXSHIFT*nsamp) - MXSHIFT*nsamp
pad = dblarr(MXSHIFT*nsamp)

FOR j = 0L, nimgs-1L DO BEGIN
   IF KEYWORD_SET(HAND_SCALE) THEN BEGIN
      arrflux[*, j] = HAND_SCALE[j]*sclflux[*, j]
      arrivar[*, j] = sclivar[*, j]/HAND_SCALE[j]^2
      arrmask[*, j] = sclmask[*, j]
      arrsky[*, j] = HAND_SCALE[j]*sclsky[*, j]
      arrnivar[*, j] = sclnivar[*, j]/HAND_SCALE[j]^2
   ENDIF ELSE IF (avg_sn LE SN_MAX_MEDSCALE AND avg_sn GT SN_MIN_MEDSCALE) $
      OR  KEYWORD_SET(MEDSCALE) THEN BEGIN
      ;; For low SNR case just scale everything to the same median flux. 
      goodpix = WHERE(sclmask[*, j] EQ 0 AND masklam, ngood)
      IF ngood EQ 0 THEN goodpix = lindgen(nspec) ;; prevent crash
      djs_iterstat, sclflux[goodpix, j], sigrej = sigrej_int, median = med_j $
                    , invvar = (sclivar[goodpix, j] >  0.0) $
                    , sigma = sig_j
      goodpix = WHERE(sclmask[*, j] EQ 0 AND masklam AND $
                      sclflux[*, j] GT 0.5*med_j, ngood)
      IF ngood EQ 0 THEN goodpix = lindgen(nspec) ;; prevent crash
      djs_iterstat, sclflux[goodpix, j], sigrej = sigrej_int, median = med_j $
                    , invvar = (sclivar[goodpix, j] >  0.0) $
                    , sigma = sig_j
      ymult[*, j] = (med_ref/med_j) < 10.0D
      arrflux[*, j] = ymult[*, j]*sclflux[*, j]
      arrivar[*, j] = sclivar[*, j]/ymult[*, j]^2
      arrmask[*, j] = sclmask[*, j]
      arrsky[*, j] = ymult[*, j]*sclsky[*, j]
      arrnivar[*, j] = sclnivar[*, j]/ymult[*, j]^2
   ENDIF ELSE IF avg_sn LE SN_MIN_MEDSCALE THEN BEGIN 
      arrflux[*, j] = sclflux[*, j]
      arrivar[*, j] = sclivar[*, j]
      arrmask[*, j] = sclmask[*, j]
      arrsky[*, j]  = sclsky[*, j]
      arrnivar[*, j] = sclnivar[*, j]
   ENDIF ELSE BEGIN
      ;; Scale the data by a polynomial if there is high enough SNR. 
      ;; Order of polynomial determined by SNR. 
      if n_elements(IN_NPOLY) EQ 0 then begin
         IF avg_sn GT 25.0 THEN npoly = 5  $ ;; is this stable??
         ELSE IF avg_sn GT 8.0 THEN npoly = 3 $ 
         ELSE IF avg_sn GE 5.0 THEN npoly = 2 $
         ELSE IF avg_sn LT 5.0 THEN npoly = 1
      endif else npoly = IN_NPOLY
      ;; grow masks by three pixels to be conservative
      polymaskflu = long_grow_mask((sclmask[*, j] EQ 0), 3)
      polymaskref = long_grow_mask(refmask, 3)
      if not keyword_set(ADDERR) then ADDERR = 0.02d
      ;; Dont allow for S/N > 50 in solve_poly_ratio fitting. This
      ;; takes care of hot pixels in invvar
      gmask_flu = sclivar[*, j]*polymaskflu GT 0
      poly_ivar_flu = gmask_flu/(1.0 / (sclivar[*, j] + (1-gmask_flu)) + $
                                 adderr^2 * (abs(sclflux[*, j]))^2)
      
      gmask_ref = refivar*polymaskref GT 0
      poly_ivar_ref = gmask_ref/(1.0 / (refivar + (1-gmask_ref)) + $
                                 adderr^2 * (abs(refflux))^2)
      if not keyword_set(FMIN) then FMIN=-100.d
      if not keyword_set(FMAX) then FMAX=1d7

      ;; If this option is set, then only fit regions with S/N >
      ;; POLY_RATIO_SN during polynomial re-scaling.
      splog, 'Polynomial re-scaling only with pixels with S/N > ', poly_ratio_sn
      IF KEYWORD_SET(POLY_RATIO_SN) THEN BEGIN
         isn = WHERE(sn_weights[*, j] GE POLY_RATIO_SN, nsn, COMP = ibad, NCOMP = nbad)
         IF nsn EQ 0 THEN message, 'S/N < 0.5 everywhere, cannot use poly_ratio_sn option'
         IF nbad GT 0 THEN BEGIN 
            poly_ivar_flu[ibad] = 0.0
            poly_ivar_ref[ibad] = 0.0
         ENDIF
      ENDIF
      solve_poly_ratio, xvector, (FMAX < sclflux[*, j] > (FMIN)) $
                        , (FMAX < refflux > (FMIN)) $
                        , poly_ivar_flu, poly_ivar_ref $
                        , npoly = npoly, nback = nback $
                        , yfit = yfit, ymult = ymult1, yadd = yadd
      ;solve_poly_ratio, xvector, (FMAX < sclflux[*, j] > (FMIN)) $
      ;                  , (FMAX < refflux > (FMIN)) $
      ;                  , sclivar[*, j]*polymaskflu $
      ;                  , refivar*polymaskref, npoly = npoly, nback = nback $
      ;                  , yfit = yfit, ymult = ymult1, yadd = yadd
      ymult1 = ((ymult1 > SCL_MIN) < SCL_MAX)
      ;; JFH 09.03.2016 apply ceiling and floor to scaling to better
      ;; treat high-z QSOs
      arrflux[*, j] = ymult1*sclflux[*, j]
      arrivar[*, j] = sclivar[*, j]/ymult1^2
      arrmask[*, j] = sclmask[*, j]
      arrsky[*, j] = ymult1*sclsky[*, j]
      arrnivar[*, j] = sclnivar[*, j]/ymult1^2
      ymult[*, j] = ymult1
      ;; At high SNR, attempt to align data taking out residual flexure 
      IF NOT KEYWORD_SET(NOSHIFT) AND ((avg_sn GE SHIFT_SN) OR KEYWORD_SET(XCORR_SKY)) THEN BEGIN
         IF KEYWORD_SET(XCORR_SKY) THEN BEGIN
            arrcorr = arrsky
            refcorr = refsky
            refcorrmask = refmask
         ENDIF ELSE BEGIN
            arrcorr = arrflux
            refcorr = refflux
            refcorrmask = refmask
         ENDELSE
         this_shift = 0
         in_lam = WHERE(newloglam GE min(loglam[*, j]) AND $
                        newloglam LE max(loglam[*, j]))
         ;; grow masks to be conservative
         maskcorr1 = long_grow_mask(arrmask[in_lam, j] EQ 0, 3)
         ;; 1=good,0=bad
         smth_corr1 = ivarsmooth(arrcorr[in_lam, j], maskcorr1, 3)
         maskcorr2  = long_grow_mask(refcorrmask[in_lam], 3)
         ;;1=good,0=bad
         smth_corr2 = ivarsmooth(refcorr[in_lam], maskcorr2, 3)
         ;; only cross-correlate pixels which are not masked in both
         goodpix = WHERE(maskcorr1 AND maskcorr2  $
                         AND smth_corr1 LT  1.0d6 $
                         AND smth_corr1 GT -100.0 $
                         AND finite(smth_corr1)   $
                         AND smth_corr2 LT  1.0d6 $
                         AND smth_corr2 GT -100.0 $
                         AND finite(smth_corr2), ngood)
         corr = c2_correlate(rebin(smth_corr1[goodpix], ngood*nsamp) $
                             , rebin(smth_corr2[goodpix], ngood*nsamp) $
                             , step, /double)
         max_corr = max(corr, jmax)
         xpeak = step[jmax]/double(nsamp)
         ;;xpeak = long_find_nminima(-corr, step/double(nsamp) $
         ;;, nfind = 1, minsep = 1 $
         ;;, ypeak = ypeak, npeak = npeak $
         ;;, errcode = errcode $
         ;;, width = 1.0, /doplot $
         ;;, xplotfit = xfit, yplotfit = yfit)
         splog, 'Measured xshift=' + $
                strcompress(string(xpeak, format = '(F7.3)'), /rem)
         IF abs(xpeak) LT 0.05 THEN splog, 'Not applying such a small shift' $
         ELSE BEGIN 
            IF abs(xpeak) EQ MXSHIFT OR KEYWORD_SET(DEBUG) THEN begin
;;                print, 'long_combspec: Your spectra may not ' + $
;;                       'be aligned to within 1.5 pixels' 
;;                print, 'long_combspec:  If this plot looks ok, then rerun with a larger value of  on..' 
               splog, 'long_combspec: Your spectra are not aligned to within', MXSHIFT $
                      , ' and you hit the boundary in the cross-correlation. You should re-running' $
                      , ' with a larger value of MXSHIFT. Type .con to continue to debug'
               stop
               print, 'This is the cross-correlatino spectrum'
               x_splot, step/double(nsamp), corr, /block $
                        , XMNX = [min(step/double(nsamp)) $
                                  , max(step/double(nsamp))] $
                        , YMNX = [min(corr), 1.0], psym1 = 10 $
                        , title = 'Cross-Corr of img # ' + $
                        strcompress(string(j), /rem) + ' with Reference.' $
                        + ' Applied xshift=' + $
                        strcompress(string(xpeak, format = '(F7.3)'), /rem) $
                        , xtwo = [xpeak], ytwo = [max_corr], psym2 = 1
               print, 'This is the shifted spectrum'
               x_splot, smth_corr1[goodpix], ytwo = smth_corr2[goodpix], /bloc
            ENDIF
               arrmask1 = interpol(arrmask[*, j], dindgen(nspec) $
                                   , dindgen(nspec) - xpeak)
               arrmask[*, j]  = (arrmask1 GT 0.5D)
               ;; interpolate over masked pixels for rebinning only
               ivar1    = djs_maskinterp(arrivar[*, j], arrmask[*, j])
               arrivar[*, j]  = interpol(ivar1,  dindgen(nspec) $
                                         , dindgen(nspec) - xpeak)* $
                                (arrmask[*, j] EQ 0) >  0.0
               arrmask[*, j] = arrmask[*, j] OR (arrivar[*, j] LE 0.0)
               flux1    = djs_maskinterp(arrflux[*, j], arrmask[*, j])
               arrflux[*, j]  = interpol(flux1, dindgen(nspec) $
                                         , dindgen(nspec) - xpeak)* $
                                (arrmask[*, j] EQ 0)
               sky1     = djs_maskinterp(arrsky[*, j], arrmask[*, j])
               arrsky[*, j]   = interpol(sky1, dindgen(nspec) $
                                         , dindgen(nspec) - xpeak)* $
                                (arrmask[*, j] EQ 0)
               nivar1   = djs_maskinterp(arrnivar[*, j], arrmask[*, j])
               arrnivar[*, j] = interpol(nivar1, dindgen(nspec) $
                                         , dindgen(nspec) - xpeak) $
                                *(arrmask[*, j] EQ 0)
               this_shift = 1
             ENDELSE
         ENDIF
   ENDELSE
    IF keyword_set(CHECK) THEN BEGIN
        IF KEYWORD_SET(THIS_SHIFT) THEN BEGIN
            x_splot, step/double(nsamp), corr, /block $
                     , XMNX = [min(step/double(nsamp)) $
                               , max(step/double(nsamp))] $
                     , YMNX = [min(corr), 1.0], psym1 = 10 $
                     , title = 'Cross-Corr of img # ' + $  ; '
                     strcompress(string(j), /rem) + ' with Reference.' $
                     + ' Applied xshift=' + $
                     strcompress(string(xpeak, format = '(F7.3)'), /rem) $
                     , xtwo = [xpeak], ytwo = [max_corr], psym2 = 1
        ENDIF
        IF KEYWORD_SET(DEBUG) THEN BEGIN
            indgood = WHERE(arrivar[*, j] GT 0.0, ngood)
            IF ngood GT 0 THEN BEGIN
                minx = min(newloglam[indgood])
                maxx = max(newloglam[indgood])
            ENDIF ELSE BEGIN
                minx = min(newloglam)
                maxx = max(newloglam)
            ENDELSE
            x_splot, 10^newloglam, refflux, ytwo = arrflux[*, j], /block $
                     , psym1 = 10, psym2 = 10 $
                     , title = 'Black = Reference, Red = Scaled' $
                     , XMNX = 10.0D^[minx, maxx] $
                     , YMNX = [ref_min, ref_max]
        ENDIF
    ENDIF
ENDFOR
outmask = round(arrmask)      ;; this converts the dbl arrmask to long

;; Only reject outliers if y
IF NOT KEYWORD_SET(NOREJ) THEN BEGIN
;avg_sn GT 3.0 AND 
    sigrej_final = sigrej
;    IF avg_sn LT 6.0 THEN sigrej_final = 5.0D $
;    ELSE sigrej_final = sigrej
    ;; Iterative rejection of outlier pixels
    niter = 5
    FOR ii = 0L, niter-1L DO BEGIN
        ;; Compute a stacked spectrum using the current mask
        imweights   = sn2_weights*(outmask EQ 0)
        newmask_now = total(outmask, 2) NE nimgs ; 1=good, 0=bad
        wght_sum    = total(imweights, 2)
        newflux_now = total(imweights*arrflux, 2)/(wght_sum + $
                                                   (wght_sum EQ 0.0))
        var =  double(arrivar NE 0.0)/(arrivar + (arrivar EQ 0.0))
        newvar =  total(imweights^2*var, 2)/(wght_sum + (wght_sum EQ 0.0))^2
        newflux_now = djs_maskinterp(newflux_now, (newmask_now EQ 0))
        newvar = djs_maskinterp(newvar, (newmask_now EQ 0))
        FOR j = 0L, nimgs-1L DO BEGIN
            ;; update noise model for rejection.  Variance = jth variance 
            ;; plus stacked spectrum variance
            var_tot = newvar + $
              double(arrivar[*, j] GT 0.0)/(arrivar[*, j] + $
                                            (arrivar[*, j] EQ 0.0))
            ivar_real = double(var_tot GT 0.0)/(var_tot + (var_tot EQ 0.0))
            ;; smooth out possible outliers in noise
            var_med = djs_median(var_tot, width = 4, boundary = 'reflect')
            var_smooth = djs_median(var_tot, width = 100, boundary = 'reflect')
            ;; conservatively always take the largest variance
            var_final = (var_med >  var_smooth) 
            ivar_final =  double(var_final GT 0.0)/$
              (var_final + (var_final EQ 0.0))
            ;; Cap S/N ratio at SN_MAX to prevent overly aggressive rejection
            SN_MAX = 20.0D
            ivar_cap = ivar_final < $
              (SN_MAX/newflux_now + (newflux_now LE 0.0))^2
            ;; adjust rejection to reflect the statistics of the distribtuion
            ;; of errors. This fixes cases where for not totally understood
            ;; reasons the noise model is not quite right and 
            ;; many pixels are rejected. 
            diff1 = arrflux[*, j]-newflux_now
            idum = where(arrmask[*, j] EQ 0, nnotmask)
            nmed_diff = nnotmask/20 > 10L
            ;; take out the smoothly varying piece 
            ;; JXP -- This isnt going to work well if the data has a bunch of
            ;; null values in it
            diff_sm = smooth(djs_median(diff1*(arrmask[*, j] EQ 0), width = nmed_diff $
                                        , boundary = 'reflect'), 5)
            chi2  = (diff1-diff_sm)^2*ivar_real
            goodchi = where(outmask[*, j] EQ 0 AND ivar_real GT 0.0 $
                            AND chi2 LE 36.0D AND masklam, ngd)
            IF ngd EQ 0 THEN goodchi = lindgen(nspec)
            ;; Is the model offset relative to the data? If so take it out
            IF NOT KEYWORD_SET(NO_OFFSET) THEN djs_iterstat, (arrflux[goodchi, j]-newflux_now[goodchi]) $
               , invvar = ivar_real[goodchi], mean = offset_mean $
               , median = offset $
            ELSE offset = 0.0
            chi2  = (arrflux[*, j]-newflux_now - offset)^2*ivar_real
            goodchi = where(outmask[*, j] EQ 0 AND ivar_real GT 0.0 $
                            AND chi2 LE 36.0D AND masklam, ngd)
            IF ngd EQ 0 THEN goodchi = lindgen(nspec)
            ;; evalute statistics of chi2 for good pixels and excluding 
            ;; extreme 6-sigma outliers
            chi2_good = chi2[goodchi]
            chi2_srt = chi2_good[sort(chi2_good)]
            ;; evaluate at 1-sigma and then scale
            gauss_prob = 1.0D - 2.0D*gaussint(-double(1.0d))
            sigind = round(gauss_prob*double(ngd))
;            sigind = 0 > (round(gauss_prob*double(ngd)) < (ngd-1))  ;; JXP 11DEC2009
            chi2_sigrej = chi2_srt[sigind]
            one_sigma = (sqrt(chi2_sigrej) > 1.0) <  5.0D
            sigrej_eff = sigrej_final*one_sigma
            chi2_cap = (arrflux[*, j]-newflux_now - offset)^2*ivar_cap
            outmask[*, j] = (arrmask[*, j] EQ 1) OR (chi2_cap GT sigrej_eff^2)
;            arr_in = arrflux[*, j] 
;            new_in = newflux_now +offset
            ;qdone = djs_reject(arr_in, new_in $
            ;                   , invvar = ivar_cap $
            ;                   , upper = sigrej_final*one_sigma $
            ;                   , lower = sigrej_final*one_sigma $
            ;                   , inmask = (arrmask[*, j] EQ 0) $
            ;                   , outmask = outmask_temp)
            ;; update the mask
;;            outmask[*, j] = (outmask_temp EQ 0) 
            ;; djs_reject returns good=1 bad=0
            IF KEYWORD_SET(CHECK) AND (ii EQ niter-1L) THEN BEGIN
                splog, 'Measured effective rejection from distribution of chi^2'
                splog, 'Instead of rejecting sigrej= ' + $
                       strcompress(string(sigrej_final), /rem) + $
                       '. Use threshold sigrej_eff= ' + $
                       strcompress(string(sigrej_eff $
                                          , format = '(F7.2)'), /rem)
                ;;gdtmp = where(outmask[*, j] EQ 0 AND ivar_real GT 0.0)
                gdtmp = where(outmask[*, j] EQ 0 AND ivar_real GT 0.0 AND masklam, ngd) ;; JFH testing
                chi = (arrflux[gdtmp, j]-newflux_now[gdtmp] - offset)* $
                  sqrt(ivar_real[gdtmp])
                resid_gauss_plot, chi, one_sigma
                ;stop
                ind_good = WHERE(arrmask[*, j] EQ 0, ngood)
                ind_bad  = WHERE(arrmask[*, j] EQ 0 AND outmask[*, j] EQ 1 $
                                 , nbad)
                mask_frac = double(nbad)/double(ngood)
                splog, 'For img ', strcompress(string(j, format = '(I3)') $
                                               , /rem) $
                       , ' nrej = ', strcompress(string(nbad, format = '(I4)') $
                                                 , /rem) $
                       , ' pixels rejected.' 
                print, '                  frac = ' $
                       , strcompress(string(mask_frac, format = '(F5.3)') $
                                     , /rem) $
                       , ' of original good pixels masked'
                indgood = WHERE(arrivar[*, j] GT 0.0, ngood)
                IF ngood GT 0 THEN BEGIN
                    minx = min(newloglam[indgood])
                    maxx = max(newloglam[indgood])
                ENDIF ELSE BEGIN
                    minx = min(newloglam)
                    maxx = max(newloglam)
                ENDELSE
                arr_in = arrflux[*, j]
                new_in = newflux_now +offset
                ;;arr_in = convol(arrflux[*, j], kernel, /EDGE_TRUN)
                ;;new_in = convol(newflux_now +offset, kernel, /EDGE_TRUNC)
                IF nbad EQ 0 THEN x_splot, 10^newloglam, new_in $
                  , ytwo = arr_in, psym1 = 10, psym2 = 10, /block $ 
                  , title = 'Black=Coadd Model, Red=' + strcompress(string(j) $
                                                                    , /rem) $
                  + 'th exposure. NO REJECTIONS', XMNX = 10.0D^[minx, maxx] $
                  , YMNX = [ref_min, ref_max] $
                  , yfou = sqrt(var_tot > 0.0), psym4 = 10 $
                ELSE x_splot, 10.0D^newloglam, new_in $
                  , ytwo = arr_in $
                  , xthr = 10.0D^newloglam[ind_bad] $
                  , ythr = arr_in[ind_bad] $
                  , yfou = sqrt(var_tot > 0.0), psym4 = 10 $
                  , psym1 = 10, psym2 = 10, psym3 = 6, /block $
                  ,  title = 'Black=Coadd Model, Red=' +  $
                              strcompress(string(j), /rem) $
                  + 'th exposure, Blue=' + strcompress(string(nbad), /rem) $  
                              + ' Rejected pixels, frac = ' + $
                              strcompress(string(mask_frac, format = '(F5.3)'), /rem) +  $
                              ' of original good pixels masked' $
                              , XMNX = 10.0D^[minx, maxx]  $
                              , YMNX = [ref_min, ref_max] 
            ENDIF
        ENDFOR
    ENDFOR
ENDIF

;; Make the final mask and also compute avsigclip combined spectra
newmask = total(outmask, 2) NE nimgs ; bad everywhere
sig2 = 1.0/(arrivar + (arrivar LE 0.0))
nused = total(outmask EQ 0, 2)
newsig2 = total(sig2*(outmask EQ 0), 2)/(nused^2 + (nused EQ 0))
newivar_sig = newmask/(newsig2 + (newsig2 LE 0.0))
;; Combine the spectra by taking a SNR weighted average. 
imweights = sn2_weights*(outmask EQ 0)
wght_sum = total(imweights, 2)
newflux = total(imweights*arrflux, 2)/(wght_sum + (wght_sum EQ 0.0))
newsky  = total(imweights*arrsky, 2)/(wght_sum + (wght_sum EQ 0.0))

IF KEYWORD_SET(CHECK) THEN BEGIN
    loadct, 0
    color_vec = lonarr(10)
    color_vec[0] = djs_icolor('white')
    color_vec[1] = djs_icolor('red')
    color_vec[2] = djs_icolor('green')
    color_vec[3] = djs_icolor('blue')
    color_vec[4] = djs_icolor('cyanxs')
    color_vec[5] = djs_icolor('magenta')
    color_vec[6] = djs_icolor('yellow')
    color_vec[7] = djs_icolor('orange')
    color_vec[8] = djs_icolor('purple')
    color_vec[9] = djs_icolor('tan')
    indgood = WHERE(total(imweights, 2) GT 0.0, ngood)
    IF ngood GT 0 THEN BEGIN
        minx = min(newloglam[indgood])
        maxx = max(newloglam[indgood])
    ENDIF ELSE BEGIN
        minx = min(newloglam)
        maxx = max(newloglam)
    ENDELSE
    plot, 10.0d^newloglam, imweights[*, 0] $
          , xrange = 10.0d^[minx, maxx] $
          , yrange = [-0.5, 1.1*max(imweights)] $
          , thick = 3, /xstyle, /ystyle, color = color_vec[0] $
          , xtitle = 'lambda', ytitle = '(S/N)^2 weight' $
          , title = 'Relative Weights Used' $
          , background = djs_icolor('black')
    FOR j = 1L, nimgs-1L DO oplot, 10.0d^newloglam, imweights[*, j] $
          , thick = 3, color = color_vec[j MOD 10]
ENDIF
;; Compute inverse variance using these weights
var =  double(arrivar NE 0.0)/(arrivar + (arrivar EQ 0.0))
newvar =  total(imweights^2*var, 2)/(wght_sum + (wght_sum EQ 0.0))^2
newivar = 1.0D/(newvar + (newvar EQ 0.0))
badflux = where(finite(newflux) NE 1 OR abs(newflux) GT 1.0d20, nbadf)
IF nbadf GT 0 THEN BEGIN
    newmask[badflux] = 0
    newflux[badflux] = 0.0
ENDIF
badivar =  where(finite(newivar) NE 1 OR abs(newivar) GT 1.0d20, nbadi)
IF nbadi GT 0 THEN BEGIN
    newmask[badivar] = 0
    newivar[badivar] = 0.0
ENDIF
newflux = newflux*newmask
newivar = newivar*newmask
;; Inverse variance without object counting photon noise. 
novar =  double(arrnivar NE 0.0)/(arrnivar + (arrnivar EQ 0.0))
newnovar = total(imweights^2*novar, 2)/(wght_sum + (wght_sum EQ 0.0))^2
newnivar = 1.0D/(newnovar + (newnovar EQ 0.0))
newnivar = newnivar*newmask

RETURN
END
