;+
; NAME:
;   long_fluxcal
;
; PURPOSE:
;   Given an extracted spectrum (in units of counts/pixel) and the
;   response function of the spectrograph, flux calibrate the spectrum
;
; CALLING SEQUENCE:
;   long_fluxcal
;
; INPUTS:
;   final_struct    - object structure containing the uncalibrated extracted
;                     spectra
;   sensfunc        - structure containing the sensitivity function
;   exptime         - exposure time for the spectra (in s)
;   airmass         - mean airmass over the course of the observations
;
; OPTIONAL INPUTS:
;  max_extrap       - Maximum delta wavelength to extrapolate the
;                     sensitivity function. [Default = 0. (Ang)]
;  /sky             - Use sky flux for the flux!  Only for MK sky
;                     modeling
;  /noextinct      - Do not apply Extinction correction
;
; OUTPUTS:
;   outfil     - Filename to contain the fluxed spectrum
;
; OPTIONAL OUTPUTS:
;  BOX_SIZE -- Size of the boxcar in native pixels (required for Sky
;              model)
;  SCIHDR
;
;
; COMMENTS:
;                   If the wavelength range of the observed spectrum
;                   differs from the wavelength range of the
;                   sensitivity function (which is likely), then the
;                   output calibrated spectrum will cover a smaller
;                   wavelength range than the input spectrum
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;                  create_struct (Goddard)
;                  struct_addtage (idlutils)
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   03-Sep-2005  Written by J. Hennawi UC Berkeley
;------------------------------------------------------------------------------
PRO long_fluxcal, scifile, sensfuncfile = sensfuncfile1 $
                  , OUTFIL = OUTFIL, SCALEFILE = SCALEFILE $
                  , WAVE = WAVE, FLUX = FLAM, SIG = FLAM_SIG $
                  , MAG = MAG, FILTER = FILTER, SCIHDR=scihdr $
                  , IN_NORM = IN_NORM, ARCHIVE = ARCHIVE $
                  , A_NORM = A_NORM, NOEXTINCT = noextinct, FRM_SCI=FRM_SCI $
                  , MAX_EXTRAP = max_extrap, SKY=sky, BOX_SIZE=box_size
                  

IF NOT keyword_set(MAX_EXTRAP) then max_extrap = 0.

IF NOT keyword_set(FRM_SCI) then begin
    scihdr = headfits(scifile)
    flux = mrdfits(scifile, 0)
    sig  = mrdfits(scifile, 1)
    wave = mrdfits(scifile, 2)
    nsig = mrdfits(scifile, 3)
;;flux = x_readspec(scifile, inflg = 2, head = scihdr, wav = wave, sig = sig)
endif else begin
    fin_strct = xmrdfits(scifile,5,finhdr, /sile)
    scihdr = xheadfits(scifile)
    if not keyword_set(SKY) then begin 
       ;; Standard approach
       flux = fin_strct[FRM_SCI-1].flux_opt
       wave = fin_strct[FRM_SCI-1].wave_opt
       sig = 1.0D/sqrt(fin_strct[FRM_SCI-1].ivar_opt)
       nsig = 1.0D/sqrt(fin_strct[FRM_SCI-1].NIVAR_OPT)
    endif else begin ;; For sky model
       if FRM_SCI LT 0 then begin
          ;; Assuming long slit
          mn = min(abs(fin_strct.xfracpos-0.46), NEW_FRM_SCI) ;; Slitb
          if mn GT 0.05 then print, 'long_fluxcal: BIG OFFSET!', mn
       endif else NEW_FRM_SCI=FRM_SCI-1L
       flux = fin_strct[NEW_FRM_SCI].sky_box
       wave = fin_strct[NEW_FRM_SCI].wave_box
       sig = 1.0D/sqrt(fin_strct[NEW_FRM_SCI].ivar_opt)
       nsig = 1.0D/sqrt(fin_strct[NEW_FRM_SCI].NIVAR_OPT)

       box_size = 2*fin_strct[NEW_FRM_SCI].box_rad* $
                  fin_strct[NEW_FRM_SCI].binning[0] ;; Native pixels of full boxcar
    endelse

    sxaddpar, scihdr, 'BITPIX', -32
    sxaddpar, scihdr, 'NAXIS', 1
    sxaddpar, scihdr, 'NAXIS1', n_elements(flux)
    sxdelpar, scihdr, 'NAXIS2'
    sxdelpar, scihdr, 'BZERO'
    sxdelpar, scihdr, 'BSCALE'
endelse
npix = n_elements(flux)


;; if no sensfunc file is passed use archived sensitivity function
IF NOT KEYWORD_SET(SENSFUNCFILE1) AND KEYWORD_SET(ARCHIVE) THEN BEGIN
    telescope = strcompress(sxpar(scihdr[*, 0], 'TELESCOP'), /rem)
    instrument = strcompress(sxpar(scihdr[*, 0], 'INSTRUME'), /rem)
    detector   =  strcompress(sxpar(scihdr[*, 0], 'DETECTOR'), /rem)
    IF strmatch(instrument, 'LRIS*') THEN BEGIN
        grism = strtrim(sxpar(scihdr[*, 0], 'GRISNAME'))
        CASE grism OF 
            '1200/3400': sensfuncfile =  GETENV('LONGSLIT_DIR') + $
              '/calib/flux/LRIS/lris_blue_1200_sensfunc.fits'
            ELSE: message, 'No archive for this setup'
        ENDCASE
    ENDIF ELSE IF strmatch(telescope, '*Gemini*') THEN BEGIN
        grating = strcompress(sxpar(scihdr[*, 0], 'GRATING'), /rem)
        CASE grating OF 
            'B1200+_G5301': sensfuncfile =  GETENV('LONGSLIT_DIR') + $
              '/calib/flux/GMOS/gmos_B1200_G5301_sensfunc.fits'
            'R150+_G5306': sensfuncfile = GETENV('LONGSLIT_DIR') + $
             '/calib/flux/GMOS/gmos_R150_G5306_sensfunc_Wolf1346.fits'
            'R400+_G5305': sensfuncfile = GETENV('LONGSLIT_DIR') + $
             '/calib/flux/GMOS/gmos_R400_GG455_sensfunc_Wolf1346.fits'
            ELSE: message, 'No archive for this setup'
        ENDCASE
    ENDIF ELSE message, 'Unrecognized detector'
ENDIF ELSE sensfuncfile = sensfuncfile1
        
mag_set = xmrdfits(sensfuncfile, 1, senshdr, /silen)
wave_min = mag_set.WAVE_MIN
wave_max = mag_set.WAVE_MAX

;; parse headers and read in extinction file
ext = long_extinct(wave, scihdr, AIRMASS = AIRMASS, EXPTIME = EXPTIME)
if keyword_set(NOEXTINCT) then ext = replicate(1./EXPTIME, n_elements(wave))
;ext = 10D^(0.4D*mag_ext*airmass)/exptime


;; Now allows for some extrapolation [JXP; May 23, 2011]
inds = WHERE( (wave GE (wave_min-MAX_EXTRAP)) AND (wave LE (wave_max + MAX_EXTRAP)))
mag_func = bspline_valu(wave[inds], mag_set)
sens = 10.0^(0.4D*mag_func)
scale = fltarr(npix)
;; KHRR 12.11.10
;; changed the next line to the line below
;scale[inds] = ext*sens
scale[inds] = sens*ext[inds]
;; scale = (extinction/exptime)*sensitivity, i.e. it is exactly what
;; is actually multiplied into the data . sens is the sensitivity func


flam = fltarr(npix)
flam_sig = fltarr(npix)
flam_nsig = fltarr(npix)
flam = flux*scale
flam_sig = sig*scale
flam_nsig = nsig*scale

; force the flux calibrated spectrum to have a given magnitude
IF KEYWORD_SET(MAG) AND NOT KEYWORD_SET(IN_NORM) THEN BEGIN
    IF NOT KEYWORD_SET(MED_WIDTH) THEN MED_WIDTH = 10L
; f_lambda= 1.0e-17 erg/cm^2/s/A = 1.0e-17*1.0e8 erg/cm^2/s/cm
; 1 Jansky = 1e-23 erg/cm^2/s/Hz
    c = 2.9979246d10
    angstrom = double(1.0d-8)
    fnu_scale =  9.1865627d-17 
;;  fnu_scale = 1.0d-8*1.0d-17/c/(3631*Jansky)
; 1.0e-17*1.0e8*(angstrom)^2/c/Jansky/3631
    LN10 = alog(10.0D)
    f_lam_med = djs_median(flam[inds], width = MED_WIDTH, boundary = 'reflect')
    f_nu = reverse((wave[inds])^2*f_lam_med)
    lognu = reverse(alog10(c/(wave[inds]*angstrom)))
    F_m = sdss_filter_nu(lognu, filter)
    mag_int = int_tabulated(lognu, LN10*f_nu*F_m)
    logA_norm = (-0.4D*mag) - alog10(fnu_scale) - alog10(mag_int)
    A_norm = 10.0D^logA_norm
    flam = A_norm*flam
    flam_sig = A_norm*flam_sig
    flam_nsig = A_norm*flam_nsig
    scale = A_norm*scale
    sens = A_norm*sens
 ENDIF ELSE IF KEYWORD_SET(IN_NORM) THEN BEGIN
    flam = in_norm*flam
    flam_sig = in_norm*flam_sig
    flam_nsig = in_norm*flam_nsig
    scale = IN_NORM*scale
    sens = IN_NORM*SENS
ENDIF
IF keyword_set(OUTFIL) THEN BEGIN
    mwrfits, flam, outfil, scihdr, /create
    mwrfits, flam_sig, outfil
    mwrfits, wave, outfil
    mwrfits, flam_nsig, outfil
    mwrfits, scale, outfil
    print, 'long_fluxcal: Final file is ', outfil
ENDIF

;; Scale is exactly what the quasar was multiplied by 
;; i.e. (includes extinction and exposure time)
;; Sens is the sensitivity function itself

RETURN
END
