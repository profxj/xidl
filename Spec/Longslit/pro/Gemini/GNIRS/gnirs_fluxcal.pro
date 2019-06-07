PRO GNIRS_FLUXCAL, scifiles, sensfuncfiles1, outfile = outfile $
                   , loglam = newloglam, flux = flux, ivar = ivar $
                   , mask = mask, check = check, box = box $
                   , HAND_SCALE = HAND_SCALE, sigrej=sigrej1

  
IF KEYWORD_SET(SIGREJ1) THEN sigrej = sigrej1 ELSE SIGREJ = 3.0D
IF NOT KEYWORD_SET(INDLAM) THEN INDLAM = 0
nfiles = n_elements(scifiles)
nimgs = 2*nfiles
exptime = dblarr(nimgs)

FOR j = 0L, nfiles-1L DO BEGIN
   obj = xmrdfits(scifiles[j], 4)
   ;; This kludge deals with outdated reductions for which the
   ;; relevant structure is in extension 5
   dims = size(obj, /dim)   
   if n_elements(dims) GT 1 THEN obj = xmrdfits(scifiles[j], 5)
   ;; Modified to 5 by JXP -- 23 July 2014
   ;; This appears only to be valid for older outdated reductions. Changed
   ;; back to 4
    hdr1 = xheadfits(scifiles[j])
    exptime[2*j:2*j+1] = double(sxpar(hdr1, 'EXPTIME'))
    IF j EQ 0 THEN BEGIN
        norders   = n_elements(obj)/2L
        nspec     = n_elements(obj[0].WAVE_OPT)
        inloglam  = dblarr(nspec, norders, nimgs)
        influx    = dblarr(nspec, norders, nimgs)
        inivar    = dblarr(nspec, norders, nimgs)
        inmask    = dblarr(nspec, norders, nimgs)
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
    ENDFOR
ENDFOR

; Take this out after fixing gnirs_reduce
inloglam = reverse(inloglam, 1)
influx = reverse(influx, 1)
inivar = reverse(inivar, 1)
scihdr = xheadfits(scifiles[0])

dims = size(influx, /dim)
nspec = dims[0]
norders = dims[1]
nimgs   = dims[2]
mask      = dblarr(nspec, norders, nimgs)

; Flux calibrate/remove telluric features spectra 
IF KEYWORD_SET(SENSFUNCFILES1) THEN BEGIN
    IF n_elements(sensfuncfiles1) EQ 1 THEN $
      sensfuncfiles = replicate(sensfuncfiles1, nimgs) $
    ELSE IF n_elements(sensfuncfiles1) EQ nfiles THEN BEGIN
        sensfuncfiles = replicate(sensfuncfiles1[0], 2)
        FOR k = 1L, nfiles-1L DO $
          sensfuncfiles = [sensfuncfiles, replicate(sensfuncfiles1[k], 2)] 
    ENDIF ELSE message, 'sensfuncfiles must be = 1 or nfiles'
    ;; Read in sensitivity function and interpolate onto new grid
    FOR j = 0L, nimgs-1L DO BEGIN
        magfunc  = xmrdfits(sensfuncfiles[j], 0)
        loglams = xmrdfits(sensfuncfiles[j], 1)
        FOR k = 0L, norders-1L DO BEGIN
            magfunc1 = interpol(magfunc[*, k], loglams, inloglam[*, k, j])
            sensfunc = 10.0D^(0.4D*magfunc1)
            scale =  sensfunc/exptime[j]
            influx[*, k, j] = influx[*, k, j]*scale
            inivar[*, k, j] = inivar[*, k, j]/(scale^2 + (scale EQ 0.0))
                                ;inmask[*, k, j] = inivar[*, k, j] GT 0.0 $
            ;  AND influx[*, k, j] LT 1000.0 AND influx[*, k, j] GT -1000.0 $
            ;  AND magfunc1 GT -10.0 AND magfunc1 LT 10.0
        ENDFOR
    ENDFOR
ENDIF
;; Define the new wavelength grid
ngrid = 5000
dloglam = 0.000127888D ; this is the average of the median dispersions
logmin = 3.777D
newloglam = logmin + dloglam*dindgen(ngrid) 
loglammask = newloglam

flux = dblarr(ngrid, norders)
ivar = dblarr(ngrid, norders)
mask = dblarr(ngrid, norders)

;; Combine the individual exposures order by order
FOR k = 0L, norders-1L DO BEGIN
    splog, 'Coadding GNIRS spectra for order # ' + $
           strcompress(string(k+ 3L), /rem)
    masklam = gnirs_ordermask(newloglam, k)
    long_combspec, reform(influx[*, k, *], nspec, nimgs) $
                   , reform(inivar[*, k, *], nspec, nimgs) $ 
                   , reform(inloglam[*, k, *], nspec, nimgs) $
                   , newloglam = newloglam $
                   , masklam = masklam, loglammask = loglammask $
                   , newflux = flux_order $
                   , newivar = ivar_order, newmask = mask_order $
                   , iref = iref, sigrej = sigrej, check = check $
                   , /NOSHIFT, MEDSCALE = MEDSCALE, HAND_SCALE = HAND_SCALE

;, MEDSCALE = (k EQ 5) 
    ;; median scale the last order
    flux[*, k] = flux_order
    ivar[*, k] = ivar_order
    mask[*, k] = mask_order
ENDFOR

IF KEYWORD_SET(OUTFILE) THEN BEGIN
    mwrfits, flux, outfile, scihdr, /create
    mwrfits, ivar, outfile
    mwrfits, mask, outfile
    mwrfits, newloglam, outfile
ENDIF

RETURN
END


