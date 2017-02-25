PRO NIRI_FLUXCAL, scifiles, sensfuncfiles1 $
                  , loglam = newloglam, flux = newflux, ivar = newivar $
                  , mask = mask,  objlam = objlam, indlam = indlam $
                  , OBJID = OBJID2, outfile = outfile $
                  , ARRFLUX = FLAM_ARR, ARRIVAR = FLIV_ARR $
                  , SIGFLUX = FLUX_SIG, SIGIVAR = IVAR_SIG, BOX = BOX $
                  , EXPTIME = EXPTIME, CHECK = CHECK, MEDSCALE = MEDSCALE  $
                  , HAND_SCALE = HAND_SCALE $
                  , LAM_MASK_MIN = LAM_MASK_MIN, LAM_MASK_MAX = LAM_MASK_MAX $
                  , IN_NPOLY = IN_NPOLY, ARR_OBJID = ARR_OBJID


IF NOT KEYWORD_SET(INDLAM) THEN INDLAM = 0
nfiles = n_elements(scifiles)

;; OBJIDVEC not specificed is the default for tellurics and we
;; then take the brightest object (or pair of objects for AB
;; sequences) on the slit
IF n_elements(objid2) EQ 0 THEN OBJID_VEC = -99 $
ELSE IF n_elements(objid2) EQ 1 THEN objid_vec = replicate(objid2-1L, nfiles) $
ELSE objid_vec = objid2 - 1L

scihdr0 = headfits(scifiles[0])
instrument = strtrim(sxpar(scihdr0, 'CURRINST'))
instrument2 = strcompress(sxpar(scihdr0, 'INSTRUME'), /rem)
;; ISAAC data in microns
;;IF strmatch(instrument2, '*ISAAC*') OR strmatch(instrument2, '*SOFI*') OR $
;;  strmatch(instrument2, '*SINFONI*') THEN wave_units = 1d4 ELSE wave_units = 1.0d

IF strmatch(instrument2, '*ISAAC*') OR $
   strmatch(instrument2, '*SINFONI*') OR strmatch(instrument2, '*SOFI*') THEN $
      wave_units = 1d4 ELSE wave_units = 1.0d

IF strmatch(instrument2, '*LUCI*') OR strmatch(instrument2, '*SOFI*') THEN BEGIN 
   dindx = 3
   IF OBJID_VEC[0] EQ -99 THEN nperfile = 2 ELSE nperfile = 1
   ntot = nfiles*nperfile
;; SINFONI apertures specified
ENDIF ELSE IF strmatch(instrument2, '*SINFONI*') AND KEYWORD_SET(ARR_OBJID) THEN BEGIN
   itot = WHERE(arr_objid GT 0, ntot)
   dindx = 4
   nperfile = -1
   IF nfiles EQ 1 THEN arr_objid = reform(arr_objid, 1, ntot)
ENDIF ELSE BEGIN 
   dindx = 4
   nperfile = 1
   ntot = nfiles*nperfile
ENDELSE

file_indx = lonarr(ntot)

ind_start = 0L
; Read in all spectra
FOR j = 0L, nfiles-1L DO BEGIN
   obj = xmrdfits(scifiles[j], dindx)
   nobj = n_elements(obj)
   IF j EQ 0 THEN BEGIN
      exptime   = dblarr(ntot)
      nspec     = n_elements(obj[0].WAVE_OPT)
      inloglam  = dblarr(nspec, ntot)
      influx    = dblarr(nspec, ntot)
      inivar    = dblarr(nspec, ntot)
;        inmask    = dblarr(nspec, nfiles)
      flux      = dblarr(nspec, ntot)
      ivar      = dblarr(nspec, ntot)
   ENDIF
   scihdr1 = headfits(scifiles[j])
   IF strmatch(instrument, 'NIRSPEC') THEN BEGIN   ;; NIRSPEC
      itime = float(strtrim(sxpar(scihdr1, 'ITIME')))
      coadds = long(strtrim(sxpar(scihdr1, 'COADDS')))
      extime_now = itime*float(coadds)
   ENDIF ELSE IF strmatch(instrument2, '*SINFONI*') OR strmatch(instrument2, '*SOFI*') $
   THEN BEGIN
      NDIT = float(esopar(scihdr1, 'HIERARCH ESO DET NDIT '))
      exptime1 = double(sxpar(scihdr1, 'EXPTIME'))
      exptime_now =  NDIT*exptime1                                
   ENDIF ELSE exptime_now = double(sxpar(scihdr1, 'EXPTIME'))    ;; NIRI
   IF OBJID_VEC[0] EQ -99 THEN BEGIN
      sn_vec = fltarr(nobj)
      FOR kk = 0L, nobj-1L DO BEGIN
         djs_iterstat, obj[kk].flux_opt*sqrt(obj[kk].ivar_opt) $
                       , mean = mean_sn, median = median_sn
         sn_vec[kk] = median_sn
      ENDFOR
      CASE nperfile OF
         1: max_obj = max(sn_vec, objid)
         2: BEGIN
            isort = reverse(sort(sn_vec))
            objid = isort[0:1]
         END
         -1: BEGIN
            ithis = where(arr_objid[j, *] GT 0, nthis)
            objid = reform(arr_objid[j, ithis]-1L, nthis)
         END
         ELSE: message, 'More than two tellurics per file not supported'
      ENDCASE
   ENDIF ELSE objid = objid_vec[j]
   CASE n_elements(objid) OF
      0: message, 'Error with objid processing'
      1: BEGIN
         exptime[j] = exptime_now
         IF KEYWORD_SET(BOX) THEN BEGIN
            loglam_temp      = alog10(obj[OBJID].WAVE_BOX*wave_units)
            inloglam[*, j]   = loglam_temp
            influx[*, j]     = obj[OBJID].FLUX_BOX
            inivar[*, j]     = obj[OBJID].IVAR_BOX
         ENDIF ELSE BEGIN
            loglam_temp      = alog10(obj[OBJID].WAVE_OPT*wave_units)
            inloglam[*, j]   = loglam_temp
            influx[*, j]     = obj[OBJID].FLUX_OPT
            inivar[*, j]     = obj[OBJID].IVAR_OPT
         ENDELSE
      END
      2: BEGIN
         exptime[2*j] = exptime_now
         exptime[2*j+1] = exptime_now
         IF KEYWORD_SET(BOX) THEN BEGIN
            loglam_temp      = alog10(obj[OBJID[0]].WAVE_BOX*wave_units)
            inloglam[*, 2*j]   = loglam_temp
            influx[*, 2*j]     = obj[OBJID[0]].FLUX_BOX
            inivar[*, 2*j]     = obj[OBJID[0]].IVAR_BOX

            loglam_temp      = alog10(obj[OBJID[1]].WAVE_BOX*wave_units)
            inloglam[*, 2*j+1]   = loglam_temp
            influx[*, 2*j+1]     = obj[OBJID[1]].FLUX_BOX
            inivar[*, 2*j+1]     = obj[OBJID[1]].IVAR_BOX
            
         ENDIF ELSE BEGIN
            loglam_temp      = alog10(obj[OBJID[0]].WAVE_OPT*wave_units)
            inloglam[*, 2*j]   = loglam_temp
            influx[*, 2*j]     = obj[OBJID[0]].FLUX_OPT
            inivar[*, 2*j]     = obj[OBJID[0]].IVAR_OPT

            loglam_temp      = alog10(obj[OBJID[1]].WAVE_OPT*wave_units)
            inloglam[*, 2*j+1]   = loglam_temp
            influx[*, 2*j+1]     = obj[OBJID[1]].FLUX_OPT
            inivar[*, 2*j+1]     = obj[OBJID[1]].IVAR_OPT
         ENDELSE
      END
      ;; This is the case of many objects
      ELSE: BEGIN
         ind_end = ind_start + nthis-1L
         exptime[ind_start:ind_end] = exptime_now
         IF KEYWORD_SET(BOX) THEN BEGIN
            loglam_temp      = alog10(obj[objid].WAVE_BOX*wave_units)
            inloglam[*, ind_start:ind_end]   = loglam_temp
            influx[*,ind_start:ind_end]     = obj[objid].FLUX_BOX
            inivar[*, ind_start:ind_end]     = obj[objid].IVAR_BOX
            file_indx[ind_start:ind_end] = j
         ENDIF ELSE BEGIN
            loglam_temp      = alog10(obj[objid].WAVE_OPT*wave_units)
            inloglam[*, ind_start:ind_end]   = loglam_temp
            influx[*, ind_start:ind_end]     = obj[objid].FLUX_OPT
            inivar[*, ind_start:ind_end]     = obj[objid].IVAR_OPT
         ENDELSE
         ind_start = ind_end + 1L
      END
   ENDCASE
ENDFOR

;; this is the output header we will use
scihdr = headfits(scifiles[indlam])

IF KEYWORD_SET(SENSFUNCFILES1) THEN BEGIN
   IF n_elements(sensfuncfiles1) EQ 1 THEN $
      sensfuncfiles = replicate(sensfuncfiles1, ntot) $
   ELSE IF n_elements(sensfuncfiles1) EQ nfiles THEN BEGIN
      CASE nperfile OF
         1: sensfuncfiles = sensfuncfiles1
         2: sensfuncfiles = reform(transpose([[sensfuncfiles], [sensfuncfiles]]), 2*nfiles)
         -1: BEGIN
            sensfuncfiles = stararr(ntot)
            sensfuncfiles = sensfuncfiles1[file_indx]
         END
         ELSE: message,  'Problem with nperfile'
      ENDCASE
   ENDIF ELSE message, 'sensfuncfiles must be = 1 or nfiles'
   ;; Read in sensitivity function and interpolate onto new grid
   dims = size(influx, /dim)
   nspec = dims[0]
;;  Flux calibrate each exposure
   FOR j = 0L, ntot-1L DO BEGIN
      magfunc1  = mrdfits(sensfuncfiles[j], 0)
      loglam_sens1   = mrdfits(sensfuncfiles[j], 1)
      magfunc = interpol(magfunc1, loglam_sens1, inloglam[*, j])
      sensfunc = 10.0D^(0.4D*magfunc)
      scale = sensfunc/exptime[j]
      influx[*, j] = influx[*, j]*scale
      inivar[*, j] = inivar[*, j]/(scale^2 + (scale EQ 0.0))
   ENDFOR
ENDIF

;; TESTING
;;in_npoly = 5
long_combspec, influx, inivar, inloglam $
               , newloglam = newloglam, newflux = newflux $
               , newivar = newivar, newmask = newmask $
               , iref = indlam, SIGREJ = SIGREJ, CHECK = CHECK $
               , /NOSHIFT, MEDSCALE = MEDSCALE, NOSHARP = NOSHARP $
               , LAM_MASK_MIN = LAM_MASK_MIN, LAM_MASK_MAX = LAM_MASK_MAX $
               , HAND_SCALE = HAND_SCALE, IN_NPOLY = IN_NPOLY

;long_combspec, flam_arr, fliv_arr, inloglam $
;               , newloglam = newloglam, newflux = flux $
;               , newivar = ivar, newmask = mask $
;               , iref = iref, SIGREJ = SIGREJ $
;               , PLT_SCALE = 1, PLT_REJ = 1, /NOSHIFT, /MEDSCALE

;; Write copmbined spectrum out to a file
newlam = 10.0D^newloglam
IF keyword_set(OUTFILE) THEN BEGIN
    sxaddpar, scihdr, 'BITPIX', -32
    sxaddpar, scihdr, 'NAXIS', 1
    sxaddpar, scihdr, 'NAXIS1', n_elements(newflux)
    sxdelpar, scihdr, 'NAXIS2'
    sxdelpar, scihdr, 'BZERO'
    sxdelpar, scihdr, 'BSCALE'
    mwrfits, newflux, outfile, scihdr, /create
    giv = where(newivar GT 0., ngiv)
    sig = 0*newivar - 1.0D
    sig[giv] = 1./sqrt(newivar[giv])
    mwrfits, sig, outfile
    mwrfits, 10.0d^newloglam, outfile
    print, 'long_coadd: Final file is ', outfile
ENDIF



; IF KEYWORD_SET(OUTFILE) THEN BEGIN
;     mwrfits, flam_arr, outfile, scihdr, /create
;     mwrfits, fliv_arr, outfile
;     mwrfits, (fliv_arr LE 0.0), outfile
;     mwrfits, loglam, outfile
;     mwrfits, flux, outfile
;     mwrfits, ivar, outfile
;     mwrfits, mask, outfile
;     mwrfits, flux_sig, outfile
;     mwrfits, ivar_sig, outfile
; ENDIF


RETURN
END


