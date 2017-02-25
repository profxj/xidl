FUNCTION LONG_COADD2D_OPT_WEIGHTS, scifiles, nbright, IREF = IREF $
                                   , YMULT = YMULT
 

  IF n_elements(iref) EQ 0 THEN iref = 0
  nimgs = n_elements(scifiles)
;; Boxcar extract the model flux image to find the 5 brightest objects
;; on the mask
  obj_model = xmrdfits(scifiles[iref], 3)
  objstruct = xmrdfits(scifiles[iref], 5)
  box_rad = objstruct[0].box_rad
  model_box = extract_boxcar(obj_model, objstruct.xpos $
                             , objstruct.ypos, radius = box_rad)
  flux_avg = djs_avsigclip(model_box, 1)
  isort = reverse(sort(flux_avg))
  ibright = isort[0:nbright-1L]
  flux_avg = flux_avg[ibright]
  bright_slitid  = objstruct[ibright].SLITID
  bright_fracpos = objstruct[ibright].XFRACPOS
  bright_objid   = objstruct[ibright].OBJID
  
  objid_arr = lonarr(nbright, nimgs)
  flux_arr = fltarr(nbright, nimgs)
  objid_arr[0:nbright-1L, iref] = bright_objid
  flux_arr[0:nbright-1L, iref] = flux_avg
;; Now loop over all other exposures and find the objects in question
  FOR iimg = 0L, nimgs-1L DO BEGIN
     IF iimg EQ IREF THEN CONTINUE
     objstruct = xmrdfits(scifiles[iimg], 5)
     obj_model = xmrdfits(scifiles[iimg], 3)
     objstruct = xmrdfits(scifiles[iimg], 5)
     model_box = extract_boxcar(obj_model, objstruct.xpos $
                                , objstruct.ypos, radius = box_rad)
     flux_avg = djs_avsigclip(model_box, 1)
     FOR iobj = 0L, nbright-1L DO BEGIN
        islit = where(objstruct.SLITID EQ bright_slitid[iobj], nslit)
        IF nslit EQ 0 THEN BEGIN
           objid_arr[iobj, iimg] = -1L
           flux_arr[iobj, iimg] = 0.0
           CONTINUE
        ENDIF ELSE BEGIN
           min_frac = min(abs(objstruct[islit].XFRACPOS-bright_fracpos[iobj]), kind)
           objid_arr[iobj, iimg] = objstruct[islit[kind]].OBJID
           flux_arr[iobj, iimg] = flux_avg[islit[kind]]
        ENDELSE
     ENDFOR
  ENDFOR
  trim = lonarr(nbright) + 1L
  FOR iobj = 0L, nbright-1L DO BEGIN
     ibad = WHERE(objid_arr[iobj, *] EQ -1L, nbad)
     IF nbad GT 0 THEN trim[iobj] = 0L
  ENDFOR
;; trim out bad objects
  igood = WHERE(trim, ngood)
  IF ngood NE nbright THEN splog, 'Of the nbright=', nbright $
                                  , ' objects considered nbad=', (nbright-ngood) $
                                  , ' could not be located'
  IF ngood EQ 0 THEN message, 'Problem, no good objects. Must be a bug'
  objid_arr = objid_arr[igood, *]
  flux_arr = flux_arr[igood, *]
  nbright = ngood
  
  flux_arr_avg = total(flux_arr, 2)/double(nimgs)
;; re-order based on average brightness
  isort = reverse(sort(flux_arr_avg))
  flux_arr = flux_arr[isort, *]
  objid_arr = objid_arr[isort, *]
  
;; Now co-add all these objects to evaluate mean SN2 for weights
  sn2_arr = fltarr(nbright, nimgs)
  FOR iobj = 0L, nbright-1L DO BEGIN
     long_coadd, scifiles, objid_arr[iobj, *], mean_sn2 = mean_sn2 $
                 , ymult = ymult1, iref = iref
     sn2_arr[iobj, *] = mean_sn2
     IF iobj EQ 0 THEN ymult = ymult1
  ENDFOR
;; Now normalize all the weights
  norm = total(sn2_arr, 2)
  weight_arr = sn2_arr/(norm # replicate(1.0, nimgs))
  weights = total(weight_arr, 1)/double(nbright)
  weights = weights/total(weights)
  
  RETURN, weights
END
