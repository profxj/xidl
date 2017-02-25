;; Based on GNIRS_FLUXCAL
;; Combines the data using long_combspec
;;    /EX_BOTH -- Both frames were extracted together (not the
;;                default)
;;  OBJ_NAM=  -- If used, then scifiles is the plan file and OBJ_NAM
;;               triggers all files with that name [preferred for
;;               scripting]
PRO TSPEC_FLUXCAL, in_scifiles, sensfuncfiles1, outfile = outfile $
                   , loglam = newloglam, flux = flux, ivar = ivar $
                   , mask = mask, check = check, box = box, EX_BOTH=ex_both $
                   , MEDSCALE = MEDSCALE, HAND_SCALE = HAND_SCALE $
                   , OBJ_NAM=obj_nam

;IF NOT KEYWORD_SET(ORDER_FIL) THEN stop
;; Generate order mask
;tset = xmrdfits(order_fil,1)
;omask = tspec_ordermask(tset, order_vec=ovec)

  if not keyword_set(OBJ_NAM) then begin
     IF NOT KEYWORD_SET(INDLAM) THEN INDLAM = 0
     scifiles = in_scifiles
  endif else begin
     ;; Read plan file
     planstr = yanny_readone(in_scifiles, hdr=planhdr, /anony) 
     ;; Parse 
     gdf = where(planstr.flavor EQ 'object' and $
                 planstr.target EQ OBJ_NAM, ngdf)
     if ngdf EQ 0 then begin
        print, 'tspec_fluxcal: No objects with name ', obj_nam
        return
     endif
     tmp_files = 'sci-'+planstr[gdf].filename
     nfil = n_elements(tmp_files)
     flg = bytarr(nfil)
     ;; Which ones actually exist??
     for ii=0L,nfil-1 do begin
        ;; Strip gz if need be
        i1 = strpos(tmp_files[ii], '.fits')
        newfil = strmid(tmp_files[ii], 0, i1+5)
        tmp_files[ii] = newfil
        ;; Does it exist
        fil = file_search(newfil+'*', count=nf)
        if nf EQ 1 then flg[ii] = 1B
     endfor
     scifiles = tmp_files[where(flg)]
  endelse
  nfiles = n_elements(scifiles)
  if keyword_set(EX_BOTH) then SIDX = 2 else SIDX = 1 
  nimgs = SIDX*nfiles 
  exptime = dblarr(nimgs)

FOR j = 0L, nfiles-1L DO BEGIN
    obj = mrdfits(scifiles[j], 4)
    hdr1 = headfits(scifiles[j])
    ;stop ;; The following needs to be fixed!
    ;exptime[2*j:2*j+1] = 600.
    if keyword_set(EX_BOTH) then begin
       exptime[2*j:2*j+1] = double(sxpar(hdr1, 'EXPTIME'))
    endif else begin
       exptime[j] = double(sxpar(hdr1, 'EXPTIME'))
    endelse
    ;; Initialize
    IF j EQ 0 THEN BEGIN
        norders   = n_elements(obj)/SIDX
        nspec     = n_elements(obj[0].WAVE_OPT)
        inloglam  = dblarr(nspec, norders, nimgs)
        influx    = dblarr(nspec, norders, nimgs)
        inivar    = dblarr(nspec, norders, nimgs)
        inmask    = dblarr(nspec, norders, nimgs)
    ENDIF
    FOR k = 0L, norders-1L DO BEGIN
        IF KEYWORD_SET(BOX) THEN BEGIN
            loglam_temp_pos       = alog10(obj[SIDX*k].WAVE_BOX)
            ;loglam_temp_neg       = alog10(obj[2*k+1L].WAVE_BOX)
            inloglam[*, k, SIDX*j]   = loglam_temp_pos
            influx[*, k, SIDX*j]     = obj[SIDX*k].FLUX_BOX
            inivar[*, k, SIDX*j]     = obj[SIDX*k].IVAR_BOX
            if keyword_set(EX_BOTH) then begin
               inloglam[*, k, 2*j+1] = loglam_temp_pos
               influx[*, k, 2*j+1]   = obj[2*k+1L].FLUX_BOX
               inivar[*, k, 2*j+1]   = obj[2*k+1L].IVAR_BOX
            endif
        ENDIF ELSE BEGIN
            loglam_temp_pos       = alog10(obj[SIDX*k].WAVE_OPT)
            ;loglam_temp_neg       = alog10(obj[2*k+1L].WAVE_OPT)
            inloglam[*, k, SIDX*j]   = loglam_temp_pos
            influx[*, k, SIDX*j]     = obj[SIDX*k].FLUX_OPT
            inivar[*, k, SIDX*j]     = obj[SIDX*k].IVAR_OPT
            if keyword_set(EX_BOTH) then begin
               inloglam[*, k, 2*j+1] = loglam_temp_pos
               influx[*, k, 2*j+1]   = obj[2*k+1L].FLUX_OPT
               inivar[*, k, 2*j+1]   = obj[2*k+1L].IVAR_OPT
            endif
        ENDELSE 
    ENDFOR
ENDFOR

scihdr = headfits(scifiles[0])

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
       stop ;; Need to deal with EX_BOTH here (JXP)
       sensfuncfiles = replicate(sensfuncfiles1[0], 2)
       FOR k = 1L, nfiles-1L DO $
          sensfuncfiles = [sensfuncfiles, replicate(sensfuncfiles1[k], 2)] 
    ENDIF ELSE message, 'sensfuncfiles must be = 1 or nfiles'
    ;; Read in sensitivity function and interpolate onto new grid
    FOR j = 0L, nimgs-1L DO BEGIN
        magfunc  = mrdfits(sensfuncfiles[j], 0)
        loglams = mrdfits(sensfuncfiles[j], 1)
        FOR k = 0L, norders-1L DO BEGIN
            magfunc1 = interpol(magfunc[*, k], loglams, inloglam[*, k, j])
            sensfunc = 10.0D^(0.4D*magfunc1)
            scale =  sensfunc/exptime[j]
            influx[*, k, j] = influx[*, k, j]*scale
            inivar[*, k, j] = inivar[*, k, j]/(scale^2 + (scale EQ 0.0))
            ;; Mask bad ivar values
            ;badiv = where(inivar[*,k,j] GT 1e10,  nbad)
            ;if nbad GT 0 then inivar[badiv,k,j] = 0.
                                ;inmask[*, k, j] = inivar[*, k, j] GT 0.0 $
            ;  AND influx[*, k, j] LT 1000.0 AND influx[*, k, j] GT -1000.0 $
            ;  AND magfunc1 GT -10.0 AND magfunc1 LT 10.0
        ENDFOR
    ENDFOR
 ENDIF

;; Define the new wavelength grid
;;   TripleSpec has a dispersion of 39.7km/s/pixel
ngrid = 8400L
dloglam = 0.000130D/alog(10) ; this is the average of the median dispersions
logmin = alog10(0.9400)  ;; microns

;;
newloglam = logmin + dloglam*dindgen(ngrid) 
loglammask = newloglam

flux = dblarr(ngrid, norders)
ivar = dblarr(ngrid, norders)
mask = dblarr(ngrid, norders)

;; Combine the individual exposures order by order
FOR k = 0L, norders-1L DO BEGIN
    splog, 'Coadding TSpec spectra for order # ' + $
           strcompress(string(7L-k), /rem)

    ;; Generate the mask
    wv0 = inloglam[*,k,0]
    gd = where(finite(wv0))
    mnwv = min(wv0, max=mxwv)
    masklam = intarr(ngrid)
    masklam[where(newloglam GE mnwv and newloglam LE mxwv)] = 1
    ;;
    long_combspec, reform(influx[*, k, *], nspec, nimgs) $
                   , reform(inivar[*, k, *], nspec, nimgs)  $ 
                   , reform(inloglam[*, k, *], nspec, nimgs) $
                   , newloglam = newloglam $
                   , masklam = masklam, loglammask = loglammask $
                   , newflux = flux_order $
                   , newivar = ivar_order, newmask = mask_order $
                   , iref = iref, sigrej = sigrej, check = check $
                   , /NOSHIFT, MEDSCALE = MEDSCALE, HAND_SCALE = HAND_SCALE

;, MEDSCALE = (k EQ 5) 
    ;; median scale the last order
    ;if k EQ 4 then stop
    flux[*, k] = flux_order
    ivar[*, k] = ivar_order
    mask[*, k] = mask_order
    ;if k EQ (norders-1) then stop
ENDFOR

IF KEYWORD_SET(OUTFILE) THEN BEGIN
   print, 'tspec_fluxcal:  Outputting fluxed spectrum ', outfile
    mwrfits, flux, outfile, scihdr, /create
    mwrfits, ivar, outfile
    mwrfits, mask, outfile
    mwrfits, newloglam, outfile
 ENDIF

RETURN
END


