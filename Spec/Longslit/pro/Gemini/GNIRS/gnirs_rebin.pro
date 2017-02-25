PRO GNIRS_PAD, inwave, influx, inivar, newloglam, dloglam, order, file $
               , outwave = outwave, outflux = outflux, outivar = outivar 

nspec = n_elements(inwave)
wv_sort = sort(inwave)
diff_srt = wv_sort-lindgen(nspec)
bad_ind = WHERE(diff_srt NE 0, nbad)
IF nbad NE 0 THEN BEGIN
    print, "WARNING: Wavelengths out of order"
    print, order, file, FORMAT = '(%"For order: %d  in file: %s")' 
    IF KEYWORD_SET(VERBOSE) THEN forprint, bad_ind, inwave[bad_ind] $
      , textout = 2
ENDIF
outpix = WHERE(inwave LT min(newloglam) OR inwave GT max(newloglam), nout)
IF nout GT 0 THEN BEGIN
    print, nout, FORMAT = '(%"DANGER: There were %d wavelengths out of range")'
    print, order, file, FORMAT = '(%"For order: %d  in file: %s")' 
    IF KEYWORD_SET(VERBOSE) THEN forprint, outpix, inwave[outpix] $
      , textout = 2
    inwave = (inwave > min(newloglam[2:*])) <  $
      max(newloglam[0:n_elements(newloglam)-2L])
ENDIF


inwave = inwave[wv_sort]
influx = influx[wv_sort]
inivar = inivar[wv_sort]

ngrid = n_elements(newloglam)
ngp_min = min(abs(min(inwave)-newloglam), ileft)
ngp_max = min(abs(max(inwave)-newloglam), iright)
nleft  = ileft-1L
nright = ngrid - nleft - nspec

inwave_left  = reverse(min(inwave) - (dindgen(nleft)+1)*dloglam)
inwave_right = max(inwave) + (dindgen(nright) + 1)*dloglam

outwave = [inwave_left, inwave, inwave_right]
outflux = [dblarr(nleft), influx, dblarr(nright)]
outivar = [dblarr(nleft), inivar, dblarr(nright)]

RETURN
END

; This routine rebins spectra onto a common wavelength grid using combine1fiber
PRO GNIRS_REBIN, scifiles, OUTFILE = OUTFILE, LOGLAM = newloglam $
                 , FLUX = flux, IVAR = IVAR, BOX = BOX

nfiles = n_elements(scifiles)
nimgs = 2*nfiles

FOR j = 0L, nfiles-1L DO BEGIN
    obj = mrdfits(scifiles[j], 5)
    IF j EQ 0 THEN BEGIN
        norders   = n_elements(obj)/2L
        nspec     = n_elements(obj[0].WAVE_OPT)
        inloglam  = dblarr(nspec, norders, 2*nfiles)
        influx    = dblarr(nspec, norders, 2*nfiles)
        inivar    = dblarr(nspec, norders, 2*nfiles)
        inmask    = dblarr(nspec, norders, 2*nfiles)
    ENDIF
    FOR k = 0L, norders-1L DO BEGIN
        IF KEYWORD_SET(BOX) THEN BEGIN
            loglam_temp_pos       = alog10(obj[2*k].WAVE_BOX)
            loglam_temp_neg       = alog10(obj[2*k+1L].WAVE_BOX)
            inloglam[*, k, 2*j]   = loglam_temp_pos
            inloglam[*, k, 2*j+1] = loglam_temp_pos
            influx[*, k, 2*j]     = obj[2*k].FLUX_BOX
            influx[*, k, 2*j+1]   = obj[2*k+1L].FLUX_BOX
            inivar[*, k, 2*j]     = obj[2*k].IVAR_BOX
            inivar[*, k, 2*j+1]   = obj[2*k+1L].IVAR_BOX
        ENDIF ELSE BEGIN
            loglam_temp_pos       = alog10(obj[2*k].WAVE_OPT)
            loglam_temp_neg       = alog10(obj[2*k+1L].WAVE_OPT)
            inloglam[*, k, 2*j]   = loglam_temp_pos
            inloglam[*, k, 2*j+1] = loglam_temp_pos
            influx[*, k, 2*j]     = obj[2*k].FLUX_OPT
            influx[*, k, 2*j+1]   = obj[2*k+1L].FLUX_OPT
            inivar[*, k, 2*j]     = obj[2*k].IVAR_OPT
            inivar[*, k, 2*j+1]   = obj[2*k+1L].IVAR_OPT
        ENDELSE 
;        inmask[*, k, 2*j]     = (inivar[*, k, 2*j] EQ 0)
;        inmask[*, k, 2*j+1]   = (inivar[*, k, 2*j+1] EQ 0)
    ENDFOR
ENDFOR

; Take this out after fixing gnirs_reduce
inloglam = reverse(inloglam, 1)
influx = reverse(influx, 1)
inivar = reverse(inivar, 1)

;logmax = 4.4025D +2.0D*dloglam
;nrebin = (logmax-logmin)/dloglam + 1 
; fixed grid wavelengths
ngrid = 4900
dloglam = 0.000127888D ; this is the average of the median dispersions
logmin = 3.777D
newloglam = logmin + dloglam*dindgen(ngrid) 

flux        = dblarr(ngrid, norders, 2*nfiles)
ivar        = dblarr(ngrid, norders, 2*nfiles)

FOR j = 0L, 2*nfiles-1L DO BEGIN
    FOR k = 0L, norders-1L DO BEGIN
;       Rebin spectra onto the newloglam grid
        inwave1 = inloglam[*, k, j]
        influx1 = influx[*, k, j]
        inivar1 = inivar[*, k, j]
        gnirs_pad, inwave1, influx1, inivar1, newloglam, dloglam $
                   , k, scifiles[j/2], outwave = outwave, outflux = outflux $
                   , outivar = outivar
        combine1fiber, outwave, outflux, newloglam = newloglam $
                       , newflux = newflux
        combine1fiber, outwave, outflux, outivar, newloglam = newloglam $
                       , newflux = newdum, newivar = newivar
        flux[*, k, j] = newflux
        ivar[*, k, j] = newivar
    ENDFOR
ENDFOR


IF KEYWORD_SET(OUTFILE) THEN BEGIN
    mwrfits, flux, outfil, /create
    mwrfits, ivar, outfile
    mwrfits, newloglam, outfile
ENDIF

RETURN
END
