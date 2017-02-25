
PRO ESI_ECHCOMBINE1D, files, outfile, IREF = IREF, NO_POLY_SCALE = NO_POLY_SCALE, CHK = CHK $
                      , LAM_MASK_MIN = LAM_MASK_MIN1, LAM_MASK_MAX = LAM_MASK_MAX1 $
                      , SHIFT_SN = SHIFT_SN1, SN_MIN_MED = SN_MIN_MED1, SN_MAX_MED1 = SN_MAX_MED1 $
                      , MXSHIFT = MXSHIFT1, NO_OFFSET = NO_OFFSET, SIGREJ = SIGREJ1 $
                      , IN_NPOLY = IN_NPOLY1, DEBUG = DEBUG $
                      , POLY_RATIO_SN = POLY_RATIO_SN1, USE_AVG_SN_WEIGHTS = USE_AVG_SN_WEIGHTS1

;IF n_elements(IREF1) GT 0 THEN IREF = IREF1 $
;ELSE IREF = 0
nspec = n_elements(files)
npix_arr = lonarr(nspec)
slit = fltarr(nspec)
;;bin_spat = lonarr(nspec)
bin_spec = lonarr(nspec)
dloglam = dblarr(nspec)
;; First sweep through all files and determine the maximum npix
FOR jj = 0L, nspec-1L DO BEGIN
   hdr = headfits(files[jj])
   binhd = strcompress(sxpar(hdr, 'BINNING'), /rem)
   binhd1 = long(strsplit(binhd, ',', /extract))
   bin_spec[jj] = binhd1[1]
   npix_arr[jj] = sxpar(hdr, 'NAXIS1')
   slitstr = sxpar(hdr, 'SLMSKNAM')
   case slitstr of 
      '0.30_arcsec': slit[jj] = 0.30
      '0.50_arcsec': slit[jj] = 0.50
      '0.75_arcsec': slit[jj] = 0.75
      '1.00_arcsec': slit[jj] = 1.00
      '1.25_arcsec': slit[jj] = 1.25
      '1.50_arcsec': slit[jj] = 1.50
      '6.00_arcsec': slit[jj] = 6.00
      else: message, 'Unrecognized slit!'
   endcase
   dloglam[jj] = double(sxpar(hdr, 'CDELT1'))
ENDFOR

;; Do we have spectra with different spectral binnings? Usually it is only 1 or
;; 2. If you binned more than 2, good luck...
i1 = WHERE(bin_spec EQ 1, nbin1)
i2 = WHERE(bin_spec EQ 2, nbin2)
IF nbin1 NE 0 AND nbin2 NE 0 THEN  BEGIN ;; Do we have hybrid binning?
   npix_arr[i1] = npix_arr[i1]/2L
   dloglam[i1] = 2.0d*dloglam[i1]
   out_spec_bin = 2L
   out_dloglam = dloglam[i2[0]] ;; Use larger wavelength grid
   hybrid_bin = 1
ENDIF ELSE IF nbin1 EQ 0 AND nbin2 NE 0 THEN BEGIN
   out_spec_bin = 2L
   out_dloglam = dloglam[0]
ENDIF ELSE IF nbin1 NE 0 AND nbin2 EQ 0 THEN BEGIN
   out_spec_bin = 1L 
   out_dloglam = dloglam[0]
ENDIF ELSE message, 'Unsupported binning'

;; Allocate arrays
nmax = max(npix_arr)
influx = fltarr(nmax, nspec)
insky = fltarr(nmax, nspec)
inivar = fltarr(nmax, nspec)
inloglam = dblarr(nmax, nspec)
FOR jj = 0L, nspec-1L DO BEGIN
   fil_sig = repstr(files[jj], '_F.fits', '_E.fits')
   fil_sky = repstr(files[jj], '_F.fits', '_S.fits')
   flux = x_readspec(files[jj], fil_sig = fil_sig, wav = wave, sig = sig)
   sky = x_readspec(fil_sky)
   loglam = alog10(wave)
   ;; Clean out NaN, bad pixels, etc. Note that fluxmax will be
   ;; problematic for bright objects!!!!
   hires_cleanspec, flux, sig, wave, fluxmin = -200.0, fluxmax = 1000.0
   IF bin_spec[jj] NE out_spec_bin THEN BEGIN
      nnow = n_elements(flux)
      nrebin = nnow/2L
      flux_rebin = fltarr(nrebin)
      loglam_rebin = dblarr(nrebin)
      sig_rebin = dblarr(nrebin)
      FOR kk = 0L, nrebin-1L DO BEGIN
         flux_temp = (flux[2*kk] + flux[2*kk+1L])/2.0D
         sig_temp = 0.5D*sqrt(sig[2*kk]^2 + sig[2*kk+1L]^2)
         IF (sig[2*kk] GT 0 AND sig[2*kk+1] GT 0) THEN BEGIN 
            flux_rebin[kk] = flux_temp
            sig_rebin[kk] = sig_temp
         ENDIF ELSE BEGIN
            flux_rebin[kk] = 0.0d
            sig_rebin[kk] = -1.0d
         ENDELSE
         loglam_rebin[kk] = (loglam[2*kk] + loglam[2*kk+1L])/2.0D
      ENDFOR
      flux = flux_rebin
      sig = sig_rebin
      loglam = loglam_rebin
   ENDIF
   ivar = (sig GT 0.0)/(sig^2 + (sig LE 0.0))
   influx[0:npix_arr[jj]-1L, jj] = flux
   insky[0:npix_arr[jj]-1L, jj] = sky
   inivar[0:npix_arr[jj]-1L, jj] = ivar
   inloglam[0:npix_arr[jj]-1L, jj] = loglam
   ;; Fill in the rest of the wavelength array with the appropriate dispersion
   IF npix_arr[jj] LT nmax THEN $
      inloglam[npix_arr[jj]:*, jj] = inloglam[npix_arr[jj]-1L, jj] + $
                                     dindgen(nmax-npix_arr[jj])*dloglam[jj] 
ENDFOR

loglam_min = min(inloglam)
loglam_max = max(inloglam)
npix_max = max(npix_arr)
;;nlam = (loglam_max - loglam_min)/(out_dloglam)
newloglam = loglam_min + dindgen(npix_max)*out_dloglam
;; Hand-scale all spectra with unit weights, i.e. disable polynomal fitting
IF KEYWORD_SET(NO_POLY_SCALE) THEN hand_scale = replicate(1.0, nspec)

long_combspec, influx, inivar, inloglam, insky = insky $
               , newloglam = newloglam $
               , innivar = innivar $
               , newflux = newflux, newivar = newivar $
               , newmask = newmask, CHECK = CHK $
               , NOSHIFT = NOSHIFT, IREF = IREF, hand_scale = hand_scale $
               , LAM_MASK_MIN = LAM_MASK_MIN1, LAM_MASK_MAX = LAM_MASK_MAX1, /XCORR_SKY $
               , SHIFT_SN = 0.0 $ ;; Always compute shifts since we use sky for ESI
               , SN_MIN_MED = SN_MIN_MED1, SN_MAX_MED = SN_MAX_MED1 $
               , MXSHIFT = MXSHIFT1, NO_OFFSET = NO_OFFSET, SIGREJ = SIGREJ1 $
               , IN_NPOLY = IN_NPOLY1, POLY_RATIO_SN = POLY_RATIO_SN1, DEBUG = DEBUG $
               , USE_AVG_SN_WEIGHTS = USE_AVG_SN_WEIGHTS1


newvar   = newmask/(newivar + (newivar EQ 0.0))
newsig = sqrt(newvar)
bad = where(newmask EQ 0, nbad)

;; Need to fix header information
IF KEYWORD_SET(OUTFILE) THEN BEGIN 
   mwrfits, newflux, outfile, /create
   mwrfits, newsig, outfile
   mwrfits, 10.0d^newloglam, outfile
ENDIF

END
