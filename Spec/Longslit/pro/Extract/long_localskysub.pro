;+
; NAME:
;   long_localskysub
;
; PURPOSE:
;   Perform sky subtraction in a given slit.  Subtraction is limited
;   to a relatively narrow region near the object.
;
; CALLING SEQUENCE:
; struct = long_localskysub( sciimg, sciivar, skyimage, piximg,
; waveimg, objstruct, thismask, skymask, objmask)
;
; INPUTS:
;
; OPTIONAL INPUTS:
; /NOLOCAL -- Turns of local sky sub (best for close pairs)
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   
; REVISION HISTORY:
;   11-Mar-2005  Written by JH + SB
;-  
;------------------------------------------------------------------------------

function long_localskysub, sciimg, sciivar, skyimage, piximg, waveimg  $
                           , ximg, objstruct, thismask, xx1, xx2, edgmask $
                           , indx = indx, bsp = bsp, skysample = skysample $
                           , niter = niter, outmask = outmask $
                           , box_rad = box_rad, objimage = objimage $
                           , modelivar = modelivar $
                           , nccd = nccd,  scihdr = scihdr $
                           , extract_struct = extract_struct $
                           , profile_struct = profile_struct $
                           , SN_GAUSS = SN_GAUSS, CHK = CHK $
                           , sigrej = sigrej1, PROF_NSIGMA=prof_nsigma1 $
                           , NOLOCAL = nolocal, STD = STD, COADD_2D = COADD_2D $
                           , VARNOOBJ = VARNOOBJ

    if NOT keyword_set(niter) then niter = 4L  
    if NOT keyword_set(bsp) THEN message, 'You must set either bsp to run long_localskysub'
    if NOT keyword_set(nccd) then nccd = 1L
    if NOT keyword_set(box_rad) then box_rad = 7
    if N_ELEMENTS(indx) EQ 0 then indx = lindgen(n_elements(objstruct))
    IF NOT KEYWORD_SET(SIGREJ1) THEN SIGREJ = 3.5D $
    ELSE sigrej = sigrej1
    IF KEYWORD_SET(PROF_NSIGMA1) THEN PROF_NSIGMA=PROF_NSIGMA1 $
    ELSE PROF_NSIGMA=fltarr(n_elements(objstruct))
   
    ;; iniatilize the outmask for the pixels on this slit
    thisind = WHERE(thismask)
    outmask[thisind] = ((thismask EQ 1) AND (sciivar GT 0))[thisind] 
    varnoobj  =  1./(modelivar + (modelivar EQ 0))
    
    nx = (size(sciimg))[1]
    ny = (size(sciimg))[2]
    xarr = findgen(nx) # replicate(1.0, ny)
    yarr = findgen(ny)## replicate(1.0, nx)
    xsize = xx2-xx1
    spatial = ximg*(xsize ## replicate(1.0, nx))
    
    nobj = n_elements(indx)
    ntot = n_elements(objstruct)
    ;; Loop over all the objects being reduced and dilate the maskwidth by 
    ;; PROF_NSIGMA for those that aplly
    FOR jj = 0L, nobj-1L DO $
       IF PROF_NSIGMA[indx[jj]] GT 0.0 THEN $
          objstruct[indx[jj]].MASKWIDTH = $
       PROF_NSIGMA[indx[jj]]*(objstruct[indx[jj]].FWHM/2.3548)
    ;; 
    o = objstruct[indx]  
    ;;------------------------------------------
    ;;  find groups of objects which may overlap
    ;;-----------------------------------------------
    i1 = 0L
    IF ARG_PRESENT(PROFILE_STRUCT) THEN BEGIN
       ;; currently I am allocating a profile image per object. 
       ;; This will become memory prohibitive if there are too many 
       ;; objects on the slit, but it is the simplest thing to do 
       ;; without dealing with the different profile image sizes
       prof_proto = {objid: 0L, slitid: 0L, mincol: 0L, maxcol:0L $
                     , MED_SN2: 0.0, profile: fltarr(nx, ny)}
       profile_struct = replicate(prof_proto, nobj)
    ENDIF

    while i1 LT nobj do begin
      group = i1
      ;; The default value of maskwidth = 3.0*FWHM = 7.05*sigma in
      ;; long_objfind, with a log(S/N) correction for bright objects
      mincols = (o[i1].xpos - o[i1].maskwidth - 1) >  xx1
      maxcols = (o[i1].xpos + o[i1].maskwidth + 1) <  xx2
      for i2=i1+1,nobj-1 do begin
         left_edge  = o[i2].xpos - o[i2].maskwidth - 1
         right_edge = o[i2].xpos + o[i2].maskwidth + 1
         touch = where(left_edge LT maxcols and o[i2].xpos GT xx1 $
                       and right_edge GT mincols, ntouch) ;; JXP -- March 2009
         if ntouch GT 0 then begin
            maxcols = (right_edge > maxcols) < xx2
            mincols = (left_edge < mincols) > xx1
;              maxcols = right_edge < xx2
              group = [group, i2]
          endif
      endfor
;;    keep for next iteration
      i1 = max(group)+1
      objwork = n_elements(group)
      scope = total(thismask EQ 1, 2)
      mincol = floor(min(mincols)) > min(where(scope))
      maxcol = ceil(max(maxcols)) < max(where(scope))
      nc = maxcol - mincol + 1L
      ipix = lindgen(nc) # replicate(1, ny) + $
        replicate(1, nc) # lindgen(ny)*nx + mincol
      mask = (thismask EQ 1) AND $
        (sciivar GT 0) AND $
        finite(sciimg) AND  $
        sciimg LT 5.0d5 AND $
        sciimg GT -1000.0d
      skymask = mask AND edgmask EQ 0 ; edge criteria
      npix = long(total(mask))
      ;; Is the the higher order npoly really necessary?? This 
      ;; makes the bspline fits take much longer and it is not clear 
      ;; it is not clear that we are gaining much here. 
      if nc GT 100 then npoly = 3L $
      else if nc GT 40 then npoly = 2L $
      else npoly = 1L
      obj_profiles = fltarr(nx*ny, objwork)
      ignoreobj = lonarr(objwork)
      sigrej_eff = sigrej
      for iiter=1L, niter do begin
        splog, 'Iteration #', iiter, ' of ', niter
        img_minsky = sciimg - skyimage
        FOR ii = 0L, objwork-1L DO BEGIN
            iobj = group[ii]
            IF iiter EQ 1L THEN BEGIN
                splog, '-------------------REDUCING--------------'
                splog, 'Finding profile for obj #' $
                       , objstruct[indx[iobj]].OBJID, ' of ', nobj
                splog, 'At  x=', djs_median(objstruct[indx[iobj]].xpos[0]) $
                       , ' on slit #', objstruct[indx[iobj]].SLITID
                splog, 'This is Object#', indx[iobj] + 1L, ' of total ', ntot
                splog, '-----------------------------------------'
                flux = extract_boxcar(img_minsky*(mask EQ 1) $
                                      , objstruct[indx[iobj]].xpos $
                                      , objstruct[indx[iobj]].ypos $
                                      , radius = box_rad)
                mvarimg = 1.0/(modelivar + (modelivar EQ 0))
                mvar_box = extract_boxcar(mvarimg*(mask EQ 1) $
                                          , objstruct[indx[iobj]].xpos $
                                          , objstruct[indx[iobj]].ypos $
                                          , radius = box_rad)
                pixtot = extract_boxcar(float(modelivar*0 + 1) $
                                        , objstruct[indx[iobj]].xpos $
                                        , objstruct[indx[iobj]].ypos $
                                        , radius = box_rad) 
                mask_box = extract_boxcar(float(mask EQ 0) $
                                          , objstruct[indx[iobj]].xpos $
                                          , objstruct[indx[iobj]].ypos $
                                          , radius = box_rad) NE pixtot
                box_denom = extract_boxcar((waveimg GT 0.0) $
                                           , objstruct[indx[iobj]].xpos $
                                           , objstruct[indx[iobj]].ypos $
                                           , radius = box_rad)
                wave = extract_boxcar(waveimg, objstruct[indx[iobj]].xpos $
                                      , objstruct[indx[iobj]].ypos $
                                      , radius = box_rad)/ $
                  (box_denom + (box_denom EQ 0))
                fluxivar = mask_box/(mvar_box + (mvar_box EQ 0))
             ENDIF ELSE BEGIN
                last_profile = reform(obj_profiles[*, ii], nx, ny)
                trace = objstruct[indx[iobj]].xpos ## replicate(1, nx)
                objmask =  ((xarr GE (trace - 2.0*box_rad)) AND $
                            (xarr LE (trace + 2.0*box_rad)))
                ;; JFH 5-20-2011 took out thismask*waveimg below
                ;; I dont think we need to mask the waveimg
                temp_str = long_extract_optimal(waveimg $
                                                , img_minsky $
                                                , sciivar*thismask $
                                                , last_profile $
                                                , outmask*objmask  $
                                                , skyimage $
                                                , objstruct[indx[iobj]].xpos $
                                                , modelivar =  $
                                                modelivar*thismask $
                                                , varnoobj =  $
                                                varnoobj*thismask $
                                                ,  box_rad = box_rad $
                                                , rn_img = rn_img)
                ;; Bad extraction?
                if temp_str.wave_opt[0] GE 0. then begin
                    ;; Good extraction
                    flux = temp_str.FLUX_OPT
                    fluxivar = temp_str.IVAR_OPT
                    wave = temp_str.WAVE_OPT
                endif
            ENDELSE
            if max(wave) gt 0 then begin
               obj_profiles[ipix, ii] = long_gprofile(img_minsky[ipix] $
                                                      , (modelivar*mask)[ipix] $
                                                      , waveimg[ipix] $         
                                                      , objstruct[indx[iobj]].xpos $
                                                      - mincol $
                                                      , wave, flux, fluxivar $
                                                      , objstruct[indx[iobj]] $
                                                      , hwidth = $
                                                      objstruct[indx[iobj]].maskwidth $
                                                      , fwhmfit = fwhmfit $
                                                      , thisfwhm = $
                                                      objstruct[indx[iobj]].fwhm $
                                                      , xnew = xnew, nccd = nccd $
                                                      , SN_GAUSS = SN_GAUSS $
                                                      , PROF_NSIGMA= $
                                                      prof_nsigma[indx[iobj]] $
                                                      , MED_SN2 = MED_SN2, /silent)
               ;; Hand apertures are currently not included in the final
               ;; object model image
               ;; JFH 09/08/2016 no longer ignore hand apertures. This
               ;; was causing a major bug and should never be
               ;; commented out. 
               ;; ignoreobj[ii] = (objstruct[indx[iobj]].HAND_AP EQ 1)
               objstruct[indx[iobj]].xpos      = xnew + mincol
               objstruct[indx[iobj]].fwhmfit   = fwhmfit
               objstruct[indx[iobj]].fwhm      = median(fwhmfit)
               mask_fact = 1.0d + 0.5d*alog10(sqrt(med_sn2 > 0.0) > 1.0)
               maskwidth = 3.0d*median(fwhmfit)*mask_fact
               IF PROF_NSIGMA[indx[iobj]] GT 0.0 THEN $
               objstruct[indx[iobj]].maskwidth = $                  
               PROF_NSIGMA[indx[iobj]]*(objstruct[indx[iobj]].FWHM/2.3548) $
               ELSE objstruct[indx[iobj]].maskwidth = maskwidth
               IF iiter EQ niter AND KEYWORD_SET(PROFILE_STRUCT) THEN BEGIN
                  ;; We are not writing these out to disk at present
                  profile_struct[iobj].OBJID = objstruct[indx[iobj]].objid
                  profile_struct[iobj].SLITID = objstruct[indx[iobj]].slitid
                  profile_struct[iobj].MINCOL = MINCOL
                  profile_struct[iobj].MAXCOL = MAXCOL
                  profile_struct[iobj].MED_SN2 = MED_SN2 
                  profile_struct[iobj].PROFILE = obj_profiles[*, ii]
                                ;IF KEYWORD_SET(profile_filename) THEN $                
                                ;  mwrfits, profile_struct, profile_filename $
                                ;  , CREATE = (indx[iobj] EQ 0), /silent
               ENDIF
            endif else begin
               print, 'long_localskysub:  Bad wavelength solution.'
               print, 'long_localskysub:  Hopefully this is an alignment box.'
               print, 'long_localskysub:  Continuing.....'
            endelse

        ENDFOR

;       Find the spatial profile of the object (at all wavelengths)
;       Use spatial profile of the sky to determine the SKY_MODEL
        sky_bmodel = 0.0
        iterbsp = 0
        WHILE total(sky_bmodel) EQ 0.0 AND iterbsp LE 5 $
          and not keyword_set(NOLOCAL) DO BEGIN
            bsp_now = (1.2^iterbsp)*bsp
            ;; if skysample is set, determine optimal break-point spacing 
            ;; directly measuring how well we are sampling of the sky. The
            ;; bsp in this case correspons to the minimum distance between 
            ;; breakpoints which we allow. 
            IF KEYWORD_SET(SKYSAMPLE) THEN BEGIN
;                sampmask = skymask
                sampmask= waveimg GT 0.0 AND thismask GT 0.0
;                FOR zz = 0L, objwork-1L DO $
;                  sampmask = sampmask $
;                  AND obj_profiles[*, zz] LT 0.05*max(obj_profiles[*, zz])
                fullbkpt = long_skybkpts(piximg, bsp_now, nx, ny $
                                         , where(sampmask))
            ENDIF ELSE BEGIN 
                ;; otherwise we use a uniform bkpt spacing of bsp
                in = where(skymask)
                pixvec = piximg[in]
                srt = sort(pixvec)
                fullbkpt = bspline_bkpts(pixvec[srt], nord = 4 $
                                         , bkspace = bsp_now, /silent) 
            ENDELSE
            ;; check to see if only a subset of the image is used. 
            ;; if so truncate since this can result in singular matrices
            xmin = min(xarr[where(thismask)])
            xmax = max(xarr[where(thismask)])
            ymin = min(yarr[where(thismask)])
            ymax = max(yarr[where(thismask)])
            isub = where(yarr GE ymin AND yarr LE ymax AND $
                         xarr GE xmin AND xarr LE xmax AND $
                         xarr GE mincol AND xarr LE maxcol)
            sortpix = sort(piximg[isub]*thismask[isub])
            ithis = where(thismask[isub])
            keep = WHERE(fullbkpt GE min(piximg[isub[ithis]]) AND $
                         fullbkpt LE max(piximg[isub[ithis]]))
            fullbkpt = fullbkpt[keep]
            sky_bmodel = long_skyoptimal(piximg[isub]*thismask[isub] $
                                         , sciimg[isub] $
                                         , (modelivar*skymask)[isub] $
                                         , obj_profiles[isub, *], sortpix $
                                         , npoly = npoly $
                                         , fullbkpt = fullbkpt $
                                         , sigrej = sigrej_eff $
                                         , obj_bmodel = obj_bmodel $
                                         , outmask = outmask1, rchi2 = rchi2 $
                                         , spatial = spatial[isub] $
                                         , ignoreobj = ignoreobj)
            iterbsp = iterbsp + 1
            IF total(sky_bmodel) EQ 0.0 AND iterbsp LE 4 THEN BEGIN 
                splog, '***************************************'
                splog, 'WARNING: bspline sky-subtraction failed'
                splog, 'Increasing bkpt spacing by 20%. Retry'
                splog, 'Old bsp = ' $
                       + strcompress(string(bsp_now $
                                            , format = '(F6.2)'), /rem) + $
                       '  New bsp = ' $
                       + strcompress(string(1.2^(iterbsp)*bsp $
                                            , format = '(F6.2)'), /rem)
                splog, '***************************************'
            ENDIF
         ENDWHILE
        rn_img = long_rdnoiseimg(nx, ny, scihdr)
        if not keyword_set(NOLOCAL) then ithispix = where(thismask[isub] EQ 1)
        IF total(sky_bmodel NE 0) THEN BEGIN
            skyimage[isub[ithispix]] = sky_bmodel[ithispix]
            objimage[isub[ithispix]] = obj_bmodel[ithispix]
            img_minsky[isub[ithispix]] = sciimg[isub[ithispix]] $
              - sky_bmodel[ithispix]
            var = abs(sky_bmodel[ithispix] + obj_bmodel[ithispix] - $
                      sqrt(2.0)*rn_img[isub[ithispix]]) $
                  + rn_img[isub[ithispix]]^2
            var_no = abs(sky_bmodel[ithispix] - $
                         sqrt(2.0)*rn_img[isub[ithispix]]) $
                     +rn_img[isub[ithispix]]^2
            ;; For weighted co-adds, the variance of the image is 
            ;; no longer equal to the image, and so the modelivar
            ;; eqn. below is not valid. However, co-adds already have 
            ;; the model noise propagated correctly in sciivar, so no
            ;; need to re-model the variance
            IF NOT KEYWORD_SET(COADD_2D) THEN BEGIN 
               modelivar[isub[ithispix]] = (var GT 0.0)/(var + (var EQ 0.0))
               varnoobj[isub[ithispix]]  = var_no
            ENDIF
            igood1 = where(skymask[isub[ithispix]], ngd1)
            ;; update the outmask for only those pixels that were fit
            ;; i.e. this prevents masking of edges in outmask
            outmask[isub[ithispix[igood1]]] = outmask1[ithispix[igood1]]
        ENDIF ELSE BEGIN
            if not keyword_set(NOLOCAL) then begin
                splog, 'ERROR: Bspline sky subtraction failed'
                IF NOT KEYWORD_SET(INBKPTS) THEN $
                  splog, '    We failed after 4 iterations of bkpt spacing. '
                splog, '       Moving on .....'
                obj_profiles = fltarr(nx*ny, objwork)
            endif
        ENDELSE
        ;;res_model = (img_minsky[ipix[ithispix]] - obj_bmodel[ithispix]) *  $
        ;;  sqrt(modelivar[ipix[ithispix]])*outmask1[ithispix]
        ;; Compute chi^2 and determine the chi^2 cut for 68% outliers
        if not keyword_set(NOLOCAL) then begin
            chi2 =  (img_minsky[isub[ithispix]] - obj_bmodel[ithispix])^2 *  $
                    modelivar[isub[ithispix]]
            ;; truncate the tail at 6-sigma if this is not a standard
            IF KEYWORD_SET(STD) THEN $
              igood = $
                WHERE(skymask[isub[ithispix]] GT 0 AND chi2 LE 100.0^2,ngd) $
            ELSE $
              igood = WHERE(skymask[isub[ithispix]] GT 0 AND chi2 LE 36.0D,ngd)
            igood2 = where(skymask[isub[ithispix]], ngd2)
            IF ngd GT 0 THEN BEGIN
               ;;DEBUG = 1
               ;;IF KEYWORD_SET(DEBUG) THEN BEGIN
               ;;   chi_model = (img_minsky[isub[ithispix]] - obj_bmodel[ithispix]) *  $
               ;;               sqrt(modelivar[isub[ithispix]])
               ;;   chi_orig = (img_minsky[isub[ithispix]] - obj_bmodel[ithispix]) *  $
               ;;              sqrt(sciivar[isub[ithispix]])
               ;;   match_hist, chi_model[igood], 100, -5.0, 5.0, /norm $
               ;;               , color = djs_icolor('black'), yr = [0.0, 0.7] $
               ;;               , xra = [-5.0, 5.0]
               ;;   match_hist, chi_orig[igood], 100, -5.0, 5.0, /norm $
               ;;               , color = djs_icolor('blue'), /oplot
               ;;   xvals = -10.0 + 0.2*findgen(101)
               ;;   ygauss1 = gauss1(xvals, [0.0, 1.0, 1.0])
               ;;   oplot, xvals, ygauss1, col = djs_icolor('red')
               ;;ENDIF
               chi2_good = chi2[igood]
                chi2_srt = chi2_good[sort(chi2_good)]
                gauss_prob = 1.0D - 2.0D*gaussint(-double(sigrej))
                sigind = round(gauss_prob*double(ngd)) < (ngd-1L)
                chi2_sigrej = chi2_srt[sigind]
                sigrej_eff = (sqrt(chi2_sigrej) >  sigrej) 
                ;; Maximum sigrej is 10 unless this is a standard. 
                IF NOT KEYWORD_SET(STD) THEN sigrej_eff = sigrej_eff < 10.0D
                splog, 'Measured effective rejection from distribution of chi^2'
                splog, 'Instead of rejecting sigrej= ' + $
                       strcompress(string(sigrej), /rem) + $
                       '. Use threshold sigrej_eff= ' + $
                       strcompress(string(sigrej_eff, format = '(F7.2)'), /rem)
                ;; Explicitly mask > sigrej outliers using the distribution of 
                ;; chi2 but only in the region that was actually fit. This
                ;; prevents e.g. excessive masking of slit edges
                outmask[isub[ithispix[igood2]]] = $
                  outmask[isub[ithispix[igood2]]] AND $
                  chi2[igood2] LT chi2_sigrej   AND $
                  sciivar[isub[ithispix[igood2]]] GT 0.0
                nrej = total(outmask[isub[ithispix[igood2]]] EQ 0)
                splog, 'Iteration ', iiter, ' rejected ', nrej, ' of ', ngd2 $
                       , ' fit pixels', format = '(a,i3,a,i8,a,i8,a)'
            ENDIF
        ENDIF
    ENDFOR

    IF KEYWORD_SET(CHK) THEN BEGIN
        qaimg = (sciimg-skyimage)*sqrt(modelivar)
        samp = lindgen(ny/4)*4
    ENDIF
    
   for ii = 0L, objwork-1L do begin
       iobj = group[ii]
       splog, 'Extracting for obj #', iobj+1, ' of ', nobj, $
              ' on slit #', objstruct[indx[iobj]].slitid $
              , ' at x=', djs_median(objstruct[indx[iobj]].xpos)
       this_profile = reform(obj_profiles[*, ii], nx, ny)
       trace = objstruct[indx[iobj]].xpos ## replicate(1, nx)
       objmask =  ((xarr GE (trace - 2.0*box_rad)) AND $
                   (xarr LE (trace + 2.0*box_rad)))
       single_struct = long_extract_optimal(waveimg $
                                            , img_minsky $
                                            , sciivar*thismask, this_profile $
                                            , outmask*objmask $
                                            , skyimage $
                                            , objstruct[indx[iobj]].xpos $
                                            , modelivar =  modelivar*thismask $
                                            , varnoobj = varnoobj*thismask $
                                            ,  box_rad = box_rad $
                                            , rn_img = rn_img)
;      single_struct.trace += mincol
       single_struct.MINCOL = MINCOL
       single_struct.MAXCOL = MAXCOL
       temp_final_struct = struct_addtags(objstruct[indx[iobj]], single_struct)
       final_struct = struct_append(final_struct, temp_final_struct)

       ;; Save to object image as necessary (i.e. NOLOCAL)
       if keyword_set(NOLOCAL) then begin
           flux_model=(replicate(1.0,nx) # single_struct.FLUX_OPT)*this_profile
           objimage = objimage + flux_model
       endif

       IF KEYWORD_SET(CHK) THEN $
         qaimg[objstruct[indx[iobj]].xpos[samp] $
               , objstruct[indx[iobj]].ypos[samp]] = -10000
   endfor

   IF KEYWORD_SET(CHK) THEN xatv, thismask*qaimg, /block $
     , min = -6.0, max = 6.0, wv = waveimg
  
endwhile

  return, final_struct
end

;; this is like boxcar extraction but does sky subtraction
;;flux_left  = extract_boxcar(img_minsky*(thismask EQ 1) $
;;                            , objstruct[indx[iobj]].xpos- $
;;                            2.*box_rad $
;;                            , objstruct[indx[iobj]].ypos $
;;                            , radius = box_rad)
;;flux_right = extract_boxcar(img_minsky*(thismask EQ 1) $
;;                            , objstruct[indx[iobj]].xpos+ $
;;                            2.*box_rad $
;;                            , objstruct[indx[iobj]].ypos $
;;                            , radius = box_rad)
;;flux = flux_center - 0.5*(flux_left + flux_right)
