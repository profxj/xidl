;+
; NAME:
;  long_coadd
;  Version 1.1
;
; PURPOSE:
;  Combine the spectra from multiple exposures taken through the same
;  mask or longslit. This routine is wrapper for the primary spectra
;  combining routine long_combspec.pro
;
; CALLING SEQUENCE:
;  LONG_COADD, infiles, objid, [OUTFIL = ]
;
; INPUTS:
;   infiles -- Names of files to coadd
;   objid   -- Object ID number for the spectrum to coadd.  This could
;              be an array if the objid changes from exposure to expsr.
;
; OPTIONAL INPUTS:
;  iref  -- Index corresponding to an exposure that defines the
;            'template' wavelength array [default: 0]
;
; OUTPUTS:
;   OUTFIL=  -- Name of the file to record coadded spectrum
;   SIGREJ=  -- Value for rejecting bad pixels [default: 2.]
;   /IRAF    -- Read wavelengths from the tag WAVE_IRAF [GMOS data]
;   EXTEN=   -- Extension in input file where the spectra structure is
;               stored [default: 5]
;   SKYFIL=  -- Name of file to record coadded sky spectrum
;
; OPTIONAL OUTPUTS:
;  FLUX=     -- Coadded flux
;  WAVE=     -- Coadded wave
;  IVAR=     -- Coadded inverse variance
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;  x_gsmooth
;
;
; REVISION HISTORY:
;   25-May-2000  Created by JH
;    4-Jun-2008  added SKYFIL option
;-
;------------------------------------------------------------------------------
;    exten=  --  Extension for the spectral file (default = 4)
;    /calib  --  Calibrated data (i.e.  use FLAM_OPT)
;    /plt_scale -- Plot the scaled data against the reference
PRO LONG_COADD, infiles, objid1, OUTFIL = OUTFIL $
                , IREF = IREF $
                , WAVE = newlam, FLUX = flux, IVAR = ivar1 $
                , ARRFLUX = flux1, ARRIVAR = IVAR $
                , IRAF = IRAF, EXTEN = exten, CALIB = calib $
                , CHECK = CHECK $
                , SIGREJ = SIGREJ1, OUTMASK = OUTMASK $
                , BOX = BOX, NOSHIFT = NOSHIFT, NOSHARP = NOSHARP $
                , NOREJ = NOREJ $
                , SKYFIL = SKYFIL, DEBUG = DEBUG,MEDSCALE=MEDSCALE $
                , HAND_SCALE = HAND_SCALE, MEAN_SN2 = MEAN_SN2, YMULT = YMULT $
                , SN_MIN_MED=SN_MIN_MED,SN_MAX_MED=SN_MAX_MAD

if  N_params() LT 2  then begin 
    print, 'Syntax - ' + $
           'long_coadd = files, objid, OUTFIL=, /BOX, /CALIB [v1.1]'
    return
endif 

nimgs = n_elements(infiles)
exptime = dblarr(nimgs)
IF NOT KEYWORD_SET(EXTEN) THEN exten = 5L 
;; JXP -- Fixed OBJID goofiness  (3/2011)
IF n_elements(objid1) EQ 1 THEN objid2 = replicate(objid1, nimgs) $
ELSE IF n_elements(objid1) EQ nimgs THEN objid2 = objid1 $
ELSE message, 'objid has wrong number of elements'

head0 = xheadfits(infiles[0])
objid = lonarr(nimgs)

FOR j = 0L, nimgs-1L DO BEGIN
    IF KEYWORD_SET(IRAF) THEN BEGIN
        obj = xmrdfits(infiles[j], 1, /sile)
        obj.WAVE_OPT = obj.WAVE_IRAF  
    ENDIF ELSE if KEYWORD_SET(BOX) then begin
        obj = xmrdfits(infiles[j], EXTEN, /sile)
        obj.WAVE_OPT = obj.WAVE_BOX
    endif else obj = xmrdfits(infiles[j], EXTEN, /sile)
    ;; Deal with objid
    a = where(obj.objid EQ objid2[j], na)
    if na NE 1 then message, 'Error in long_coadd. Could not find that objid'
    objid[j] = a[0]  ;; Ugly kludge
    ;;
    head = xheadfits(infiles[j])
    dum = long_extinct(obj[objid[j]].wave_opt, head, exptime = exptime1)
    exptime[j] = exptime1
    IF j EQ 0 THEN BEGIN
        nspec = n_elements(obj[objid[j]].WAVE_OPT)
        IF obj[objid[j]].WAVE_OPT[0L+10L] GT $
          obj[objid[j]].WAVE_OPT[nspec-1L-10L] THEN FLIP = 1
        inloglam = dblarr(nspec, nimgs)
        influx   = dblarr(nspec, nimgs)
        inivar   = dblarr(nspec, nimgs)
        innivar  = dblarr(nspec, nimgs)
        insky    = dblarr(nspec, nimgs)
        loglam = dblarr(nspec, nimgs)
    ENDIF

    ;; Save wave
;;;;;;Change 1: replace all negative wavelengths with a value of 1.0
    loglam_temp  = alog10(obj[objid[j]].wave_opt > 1.)
    inloglam[*, j] = loglam_temp

    ;; Set the flux
    if not keyword_set(CALIB) then begin
       if KEYWORD_SET(BOX) then begin
          influx[*, j]   = obj[objid[j]].FLUX_BOX
          inivar[*, j]   = obj[objid[j]].IVAR_BOX
          innivar[*, j]  = obj[objid[j]].NIVAR_BOX
          insky[*, j]    = obj[objid[j]].SKY_BOX
       endif else begin
          influx[*, j]   = obj[objid[j]].FLUX_OPT
          inivar[*, j]   = obj[objid[j]].IVAR_OPT
          innivar[*, j]  = obj[objid[j]].NIVAR_OPT
          insky[*, j]    = obj[objid[j]].SKY_OPT
       endelse
    endif else begin
       if KEYWORD_SET(BOX) then begin
          influx[*, j]   = obj[objid[j]].FLAM_BOX
          inivar[*, j]   = obj[objid[j]].FLAM_IVAR_BOX
       endif else begin
          influx[*, j]   = obj[objid[j]].FLAM_OPT
          inivar[*, j]   = obj[objid[j]].FLAM_IVAR_OPT
       endelse
    endelse
ENDFOR

IF KEYWORD_SET(FLIP) THEN BEGIN
    inloglam = reverse(inloglam)
    influx   = reverse(influx)
    inivar   = reverse(inivar)
    innivar  = reverse(innivar)
    insky    = reverse(insky)
ENDIF

long_combspec, influx, inivar, inloglam, insky = insky, innivar = innivar $
               , newloglam = newloglam, newflux = newflux $
               , newivar = newivar, newnivar = newnivar $
               , newmask = newmask, newsky = newsky $
               , iref = iref, SIGREJ = SIGREJ, CHECK = CHECK $
               , NOSHIFT = NOSHIFT, NOSHARP = NOSHARP, NOREJ = NOREJ $
               , MEDSCALE=MEDSCALE, HAND_SCALE = HAND_SCALE $
               , SN2 = SN2, YMULT = YMULT, DEBUG = DEBUG, MEAN_SN2 = MEAN_SN2 $
               , SN_MIN_MED=SN_MIN_MED,SN_MAX_MED=SN_MAX_MED

;; Im commenting out this last step, since I realize it undoes the 
;; heliocentric correction. This precludes the possibility of improving
;; the flexure correction with a higher S/N ratio sky. If I instead
;; put the heliocentric correction in last then I would not be able to 
;; to combine data from different dates. 

;; Take out any residual flexure with a higher SNR sky
; long_reduce_params, head0, [1, 1], skyfile = skyfile
; IF KEYWORD_SET(SKYFILE) AND nimgs GT 1 THEN BEGIN
;     wave = 10.0d^newloglam
;     struct = create_struct('WAVE_OPT', wave, 'WAVE_BOX', wave $
;                            , 'SKY_BOX', newsky, 'MASK_BOX', newmask  $
;                            , 'FLX_SHFT_WAV', 0.0D)
;     qafile = repstr(outfil, '.fits', '-flex.ps')
;     long_flexure, struct, skyfile, qafile = qafile
;     IF abs(struct.FLX_SHFT_WAV) LE 2.0 THEN BEGIN
;         newloglam = alog10(struct.WAVE_OPT)
;         flx_shft_wav = struct.FLX_SHFT_WAV
;     ENDIF ELSE BEGIN
;         splog, 'Coadded flexure correction too large. Not applying it'
;         FLX_SHFT_WAV = 9999.0
;     ENDELSE
; ENDIF
;; Write copmbined spectrum out to a file
newlam = 10.0D^newloglam
IF keyword_set(OUTFIL) THEN BEGIN
    sxaddpar, head0, 'NEXP', nimgs, ' number of exposures combined'
    sxaddpar, head0, 'EXPT_AVG', float(total(exptime)/double(nimgs)), $
      ' mean exposure time (s)'
    sxaddpar, head0, 'BITPIX', -32
    sxaddpar, head0, 'NAXIS', 1
    sxaddpar, head0, 'NAXIS1', n_elements(newflux)
;    IF KEYWORD_SET(SKYFILE) AND nimgs GT 1 THEN $
;      sxaddpar, head0, 'FLX_SHFT_WAV', flx_shft_wav
    sxdelpar, head0, 'NAXIS2'
    sxdelpar, head0, 'BZERO'
    sxdelpar, head0, 'BSCALE'
    mwrfits, newflux, outfil, head0, /create
    giv = where(newivar GT 0., ngiv)
    sig = 0*newivar - 1.0D
    sig[giv] = 1./sqrt(newivar[giv])
    mwrfits, sig, outfil
    mwrfits, 10.0d^newloglam, outfil
;    gniv = where(newnivar GT 0., ngniv)
    nsig = newnivar
    nsig = 0*newivar - 1.0D
    nsig[giv] = 1./sqrt(newnivar[giv])
    mwrfits, nsig, outfil
    mwrfits, sn2, outfil
    mwrfits, ymult, outfil
    print, 'long_coadd: Final file is ', outfil
ENDIF

IF KEYWORD_SET(SKYFIL) THEN BEGIN
    sxaddpar, head0, 'NEXP', nimgs, ' number of exposures combined'
    sxaddpar, head0, 'EXPT_AVG', float(total(exptime)/double(nimgs)), $
      ' mean exposure time (s)'
    sxaddpar, head0, 'BITPIX', -32
    sxaddpar, head0, 'NAXIS', 1
    sxaddpar, head0, 'NAXIS1', n_elements(newflux)
    sxdelpar, head0, 'NAXIS2'
    sxdelpar, head0, 'BZERO'
    sxdelpar, head0, 'BSCALE'
    mwrfits, newsky, skyfil, head0, /create
    sig = sqrt(newsky)          ; assum Poisson
    bd = where(newsky lt 0.,nbd)
    if nbd ne 0 then sig[bd] = 0.
    mwrfits, sig, skyfil 
    mwrfits, 10.0d^newloglam, skyfil
    mwrfits, sig, skyfil     ;just guessing at this
    print, 'long_coadd: Final file is ', skyfil
ENDIF

RETURN
END
