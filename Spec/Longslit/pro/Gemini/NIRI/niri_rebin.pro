; This routine rebins spectra onto a common wavelength grid using combine1fiber
PRO NIRI_REBIN, scifiles, loglam = loglam $
                , FLUX = flux, IVAR = IVAR $
                , NEWLOGLAM = newloglam, INDLAM = INDLAM $
                , BOX = BOX, OBJID = OBJID2, EXPTIME = EXPTIME, SCIHDR = SCIHDR

message, 'You need to fix this. We should not be rebiining with bsplines here'

nfiles = n_elements(scifiles)
IF n_elements(objid2) EQ 1 THEN objid1 = replicate(objid2, nfiles) $
ELSE objid1 = objid2
objid = objid1-1L

; Read in all spectra
FOR j = 0L, nfiles-1L DO BEGIN
    obj = mrdfits(scifiles[j], 4)
    nobj = n_elements(obj)
    IF j EQ 0 THEN BEGIN
        exptime   = dblarr(nfiles)
        nspec     = n_elements(obj[0].WAVE_OPT)
        inloglam  = dblarr(nspec, nfiles)
        influx    = dblarr(nspec, nfiles)
        inivar    = dblarr(nspec, nfiles)
;        inmask    = dblarr(nspec, nfiles)
        flux      = dblarr(nspec, nfiles)
        ivar      = dblarr(nspec, nfiles)
    ENDIF
    scihdr1 = headfits(scifiles[j])
    exptime[j] = double(sxpar(scihdr1, 'EXPTIME'))
    IF KEYWORD_SET(BOX) THEN BEGIN
        loglam_temp      = alog10(obj[OBJID[j]].WAVE_BOX)
        inloglam[*, j]   = loglam_temp
        influx[*, j]     = obj[OBJID[j]].FLUX_BOX
        inivar[*, j]     = obj[OBJID[j]].IVAR_BOX
    ENDIF ELSE BEGIN
        loglam_temp      = alog10(obj[OBJID[j]].WAVE_OPT)
        inloglam[*, j]   = loglam_temp
        influx[*, j]     = obj[OBJID[j]].FLUX_OPT
        inivar[*, j]     = obj[OBJID[j]].IVAR_OPT
    ENDELSE 
ENDFOR

IF NOT KEYWORD_SET(NEWLOGLAM) THEN BEGIN
    IF NOT KEYWORD_SET(INDLAM) THEN INDLAM = 0
    newloglam = inloglam[*, INDLAM]
    scihdr = headfits(scifiles[indlam])
ENDIF

FOR j = 0L, nfiles-1L DO BEGIN
; Rebin spectra onto the newloglam grid
    inwave1 = inloglam[*, j]
    influx1 = influx[*, j]
    inivar1 = inivar[*, j]
    combine1fiber, inwave1, influx1, newloglam = newloglam $
                   , newflux = newflux
    combine1fiber, inwave1, influx1, inivar1, newloglam = newloglam $
                   , newflux = newdum, newivar = newivar
    flux[*, j] = newflux
    ivar[*, j] = newivar
ENDFOR
loglam = newloglam

;IF KEYWORD_SET(OUTFILE) THEN BEGIN
;    mwrfits, flux, outfil, /create
;    mwrfits, ivar, outfile
;    mwrfits, newloglam, outfile
;ENDIF

RETURN
END
