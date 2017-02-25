PRO LONG_RESVEC, infiles, objid1, wave, resvec = resvec, outfil = outfil

  clight = 2.9979246e5
  if  N_params() LT 2  then begin 
     print, 'Syntax - ' + $
            'long_coadd = files, objid, OUTFIL=, /BOX, /CALIB [v1.1]'
     return
  endif 
  
  nimgs = n_elements(infiles)
  arc_fwhm_med = dblarr(nimgs)
  IF NOT KEYWORD_SET(EXTEN) THEN exten = 5L 
  IF n_elements(objid1) EQ 1 THEN objid = replicate(objid1-1L, nimgs) $
  ELSE IF n_elements(objid1) EQ nimgs THEN objid = objid1-1L $
  ELSE message, 'objid has wrong number of elements'
  
  head0 = xheadfits(infiles[0])
  FOR j = 0L, nimgs-1L DO BEGIN
     IF KEYWORD_SET(IRAF) THEN BEGIN
        obj = xmrdfits(infiles[j], 1, /sile)
        obj.WAVE_OPT = obj.WAVE_IRAF  
     ENDIF ELSE if KEYWORD_SET(BOX) then begin
        obj = xmrdfits(infiles[j], EXTEN, /sile)
        obj.WAVE_OPT = obj.WAVE_BOX
     endif else obj = xmrdfits(infiles[j], EXTEN, /sile)
     IF j EQ 0 THEN BEGIN
        nspec = n_elements(obj[objid[j]].WAVE_OPT)
        IF obj[objid[j]].WAVE_OPT[0L+10L] GT $
           obj[objid[j]].WAVE_OPT[nspec-1L-10L] THEN FLIP = 1
        inwave = dblarr(nspec, nimgs)
        arc_fwhm_fit   = dblarr(nspec, nimgs)
     ENDIF
     inwave[*, j] = obj[objid[j]].wave_opt
     arc_fwhm_fit[*, j] = obj[objid[j]].ARC_FWHM_FIT
     arc_fwhm_med[j] = obj[objid[j]].ARC_FWHM_MED
  ENDFOR

  IF KEYWORD_SET(FLIP) THEN BEGIN
     inwave = reverse(inwave)
     arc_fwhm_fit = reverse(arc_fwhm_fit)
  ENDIF

  nwave = n_elements(wave)
  resarr = fltarr(nwave, nimgs)
  resmask = lonarr(nwave, nimgs)
  mask = fltarr(nwave, nimgs)
  dwv = fltarr(nspec)
  FOR j = 0L, nimgs-1L DO BEGIN
     dwv[3:nspec-2L] = abs((inwave[*, j]-shift(inwave[*, j], 1))[2:nspec-3L]) 
     dwv[0:3] = dwv[4]
     dwv[nspec-3L:nspec-1L] = dwv[nspec-4L]
     res = djs_median(clight*(dwv*arc_fwhm_fit[*, j])/inwave[*, j] $
                      , width = 10, boundary = 'reflect')
     ind = WHERE(wave GE min(inwave[*, j]) AND wave LE max(inwave[*, j]))
     isort = sort(inwave[*, j])
     resarr[ind, j] = interpol(res[isort], inwave[isort, j], wave[ind])
     resmask[ind, j] = 1L
  ENDFOR
  IF nimgs EQ 1 THEN weight_sum = resmask $
  ELSE weight_sum = total(resmask, 2)
  IF nimgs EQ 1 THEN $
     resvec = (resarr*resmask)/(weight_sum + (weight_sum EQ 0.0)) $
  ELSE resvec = total(resarr*resmask, 2)/(weight_sum + (weight_sum EQ 0.0))
  ibad = WHERE(resvec EQ 0.0, nbad, COMPLEMENT = igood)
  IF nbad GT 0 THEN $
     resvec[ibad] = interpol(resvec[igood], wave[igood], wave[ibad])
  
IF keyword_set(OUTFIL) THEN BEGIN
    sxaddpar, head0, 'NEXP', nimgs
    sxaddpar, head0, 'BITPIX', -32
    sxaddpar, head0, 'NAXIS', 1
    sxaddpar, head0, 'NAXIS1', n_elements(resvec)
    sxdelpar, head0, 'NAXIS2'
    sxdelpar, head0, 'BZERO'
    sxdelpar, head0, 'BSCALE'
    sxaddpar, head0, 'ARC_FWHM', total(arc_fwhm_med)/double(nimgs)
    mwrfits, resvec, outfil, head0, /create
    mwrfits, resvec*0, outfil
    mwrfits, wave, outfil
    print, 'long_resvec: Resolution file is ', outfil
 ENDIF
  
RETURN
END
