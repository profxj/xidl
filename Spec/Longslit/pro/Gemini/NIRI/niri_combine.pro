PRO NIRI_COMBINE, influx1, inivar1, inmask = inmask1 $
                  , IREF = IREF1, FLUX = flux, IVAR = ivar $
                  , MASK = MASK, SIGREJ = SIGREJ1 $
                  , SIGFLUX = flux_sig, SIGIVAR = ivar_sig

influx = influx1
inivar = inivar1

SN_MIN = 5.0D
dims = size(influx, /dim)
nspec   = dims[0]
nimgs   = dims[1]

IF KEYWORD_SET(SIGREJ1) THEN sigrej = sigrej1 $
ELSE BEGIN
    if (nimgs LE 2) then sigrej = 1.0 $ 
    else if (nimgs EQ 3) then sigrej = 1.1 $
    else if (nimgs EQ 4) then sigrej = 1.3 $
    else if (nimgs EQ 5) then sigrej = 1.6 $
    else if (nimgs EQ 6) then sigrej = 1.9 $
    else sigrej = 2.0
ENDELSE


IF NOT KEYWORD_SET(INMASK1) THEN inmask = lonarr(nspec, nimgs) + 1L $
ELSE inmask = inmask1

; Measure median SNR
snr_arr  = dblarr(nimgs)
med_flux = dblarR(nimgs)
FOR j = 0L, nimgs-1L DO BEGIN
    ind_ord = WHERE(inivar[*, j] GT 0.0 AND $
                    finite(inivar[*, j]) EQ 1, nind)
    djs_iterstat, influx[ind_ord, j]*sqrt(inivar[ind_ord, j]) $
                  , sigrej = 3.0, median = median_sn
    djs_iterstat, influx[ind_ord, j], sigrej = 3.0, median = median_f
    snr_arr[j]  = median_sn
    med_flux[j] = median_f 
ENDFOR

IF NOT KEYWORD_SET(IREF1) THEN max_snr = max(snr_arr, iref) $
ELSE IREF = IREF1

finalmask = (inivar GT 0.0) $
  AND (finite(influx) EQ 1) $
  AND (finite(inivar) EQ 1) $
  AND inmask

scale    = replicate(1.0D, nspec) # (med_flux[iref]/med_flux)
fweight2 = replicate(1.0D, nspec) # (med_flux)^2
weights = inivar*fweight2

; (S/N)^2 weights. We use the median flux^2 of the whole spectrum with pixel 
; inivar weights (since noise is all sky and high SNR). 
influx = influx*scale
inivar = inivar/scale^2

;   Avsigclip the spectra to get the outmask with no weights
newflux_sig = djs_avsigclip(influx, 2, inmask = (finalmask EQ 0) $
                            , outmask = outmask, sigrej = sigrej)
newmask = total(outmask, 2) NE nimgs ; bad everywhere
sig2 = 1.0/(inivar + (inivar LE 0.0))
nused = total(outmask EQ 0, 2)
newsig2 = total(sig2*(outmask EQ 0), 2)/(nused^2 + (nused EQ 0))
newivar_sig = newmask/(newsig2 + (newsig2 LT 0.0))

;   Combine the spectra by taking a SNR weighted average. Since we are 
;   background limited the noise is known well so this is sensible
weights = weights*(outmask EQ 0)
wght_sum = total(weights, 2)
newflux = total(weights*influx, 2)/(wght_sum + (wght_sum EQ 0.0))
newivar = total(inivar, 2)
flux = newflux
ivar = newivar*newmask
mask = newmask
flux_sig = newflux_sig
ivar_sig = newivar_sig

RETURN
END
