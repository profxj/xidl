pro mage_pipe, mage, obj_id=obj_id, chk=chk $
               , clobber = clobber, sensfunc = sensfunc $
               , singles = singles, NOCR = NOCR
  
  ;; setup directories
  IF NOT FILE_TEST('Object', /DIR) THEN FILE_MKDIR, 'Object'
  IF NOT FILE_TEST('Final', /DIR)  THEN FILE_MKDIR, 'Final'
  IF NOT FILE_TEST('FSpec', /DIR)  THEN FILE_MKDIR, 'FSpec'
  ;; Pick out science targets
  isci = where(strtrim(mage.exptype) EQ 'SCIENCE' OR $
               strtrim(mage.exptype) EQ 'STD'     OR $
               strtrim(mage.exptype) EQ 'BRIGHT' AND $
               mage.obj_id GE 0, nscience)
  sciframes = mage[isci]
  obj_ids = sciframes[uniq(sciframes.obj_id)].obj_id
  ;;if we only are reducing one object
  if (keyword_set(OBJ_ID)) then begin
     obj_ids = [obj_id]
  endif

  ;; Loop over objects
  for iobj=0, n_elements(obj_ids)-1 do begin
     objframes = sciframes[where(sciframes.obj_id EQ obj_ids[iobj], nframes)]
     print, "MAGE_PIPE: Processing object ", objframes[0].object
     ;; add some header cards

     master_hdr = xheadfits(objframes[0].RAWPATH + objframes[0].FITSFILE)
     ;;master_hdr = objframes[0].hdr

     sxdelpar, master_hdr, "NAXIS1"
     sxdelpar, master_hdr, "CDELT1"
     sxdelpar, master_hdr, "BZERO"
     sxdelpar, master_hdr, "CD1_1"
     sxdelpar, master_hdr, "CRVAL1"
     sxaddpar, master_hdr, "COMMENT", "mage_pipe: v0.3, Build June 2009"
     sxaddpar, master_hdr, "COMMENT", "object reduced as: " + $
               objframes[0].exptype
     pipe_time = "mage_pipe: Reduced "+systime()
     sxaddpar, master_hdr, "COMMENT", pipe_time
     
     for iexp=0, nframes-1 do begin

        IF NOT keyword_set(CLOBBER) THEN BEGIN
           tmp = strsplit(objframes[iexp].fitsfile, 'mage', /extract)
           outfil = 'Object/ObjStr'+strtrim(tmp[0])
           
           if (file_test(outfil) EQ 1) then begin

              tt = xmrdfits(outfil, 1)
              ;reflux the reduced spectra here?
              
              IF NOT KEYWORD_SET(sensfunc) THEN BEGIN
                 print, 'Refluxing with ' + $
                        fileandpath(objframes[iexp].sensfunc)
                 mage_flux, objframes[iexp].sensfunc, tt, rej = 0.05
                 sxaddpar, master_hdr, 'COMMENT' $
                           , 'mage_pipe: Re-fluxed with ' +  $
                           strcompress(fileandpath(objframes[iexp].sensfunc) $
                                       , /rem)
              ENDIF
              
              print, "Already Reduced, adding in ", outfil
              
              if (iexp EQ 0) then begin
                 allframes = [ tt ]
              endif else begin
                 allframes = [ allframes, tt ]
              endelse
              continue
           endif
        endif
        ;; Is this a standard or bright object??
        STD = strcompress(objframes[iexp].exptype,/rem) EQ 'STD'
        BRIGHT = strcompress(objframes[iexp].exptype,/rem) EQ 'BRIGHT'
        ;; Bias subtract and flat field
        fitsfile = objframes[iexp].rawpath+objframes[iexp].fitsfile 
        pixflatfile = objframes[iexp].pixflatfile 
        illumflatfile = objframes[iexp].illumflatfile 
        
        mage_proc, fitsfile, sciimg, sciivar, pixflatfile = pixflatfile $
                   , illumflatfile = illumflatfile, hdr = scihdr

        arcfile = objframes[iexp].rawpath+objframes[iexp].arcfile
        slitfile = objframes[iexp].slitfile
        tset_slits = mrdfits(slitfile, 1)
        ;; Trace the arc + science data
        
        mage_proc, arcfile, arcimg
        piximg=mage_makescipix(arcimg,sciimg,tset_slits,chk=chk $
                               , std = std, bright = bright)
        ;; Generate a 2D skymodel
        print, "Generating the 2D sky model"
        skyimage = mage_skymodel(sciimg, sciivar, piximg=piximg $
                                 ,tset_slits=tset_slits)
        
        ;; Reject cosmic rays
        IF NOT KEYWORD_SET(NOCR) THEN BEGIN
           print, 'Cosmic ray rejection'
           sigma_psf = 3.0D/2.35482D
           ;;   Take the PSF width to be that of the spectral direction.  
           ;;   This prevents the routine from rejecting sky lines
           ;;   This is a description of the 3x3 core of the 2D PSF
           ;;   for  psf_reject_cr.pro
           ;;    
           ;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1] 
           ;;    PSFVALS[0]          1.   PSFVALS[0]
           ;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
           psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
           crmask  = psf_reject_cr(sciimg-skyimage, sciivar, psfvals $
                                   , satmask = (sciimg-skyimage) GT 8d4)
           sciivar = sciivar*(crmask EQ 0)
           ;; Do we need to sky-subtract before CR rejection???
           ;; Do we need to mask these out in the sciimg??
        ENDIF        
        ;; Generate a wavelenth solution for nearest arc
        print, "Generating wavelength solution"
        print, "Arc file: ", arcfile
        orderfile = objframes[iexp].orderfile
        ordr_str = mrdfits(orderfile, 1)

        mage_arc, arcfile, ordr_str=ordr_str, outarc=arcimg_fil, /clobber
        waveimg = xmrdfits(arcimg_fil)
        ;; Find and trace objects
        IF NOT KEYWORD_SET(STD) THEN BEGIN
           filstd = strcompress(objframes[iexp].stdfile, /rem)
           IF x_chkfil(filstd+'*') EQ 0 then begin
              print, 'mage_findobj: STD Obj file does not exist or' + $
                     ' obj_fil tag not set! Check it..', filstd
              print, "ERROR: Go back and enter your standard star file name into the GUI box (Sens func tab)"
              RETURN
           ENDIF ELSE filstd=filstd[0]
        ENDIF ELSE filstd = 0
        print, "Finding objects on ", fitsfile
        obj_strct = mage_findobj(sciimg-skyimage, sciivar, waveimg, tset_slits $
                                 , filstd = filstd, chk = chk)
        
        ;; Optimal extraction
        tmp = strsplit(objframes[iexp].fitsfile, 'mage', /extract)
        outfil = 'Final/f_'+objframes[iexp].fitsfile
        mage_echextobj, sciimg, sciivar, scihdr, skyimage $
                        , piximg, waveimg, tset_slits $
                        , obj_strct, outfil = outfil $
                        , box_rad = box_rad, STD = STD, CHK = CHK 
        
        ;; Set some items in the obj_strct
        obj_strct.img_fil = objframes[iexp].rawpath+objframes[iexp].fitsfile
        obj_strct.arc_fil = objframes[iexp].rawpath+objframes[iexp].arcfile
        obj_strct.exp=float(sxpar(scihdr,'EXPTIME'))
        
        ;; Flux calibrate spectrum
        IF NOT KEYWORD_SET(STD) THEN $
           mage_flux, objframes[iexp].sensfunc, obj_strct, rej = 0.05

        sxaddpar, master_hdr, 'COMMENT', 'mage_pipe: Fluxed with ' + $
                  strcompress(objframes[iexp].sensfunc, /rem)
        objfil = 'Object/ObjStr'+strcompress(tmp[0],/rem)
        ;;write out object file
        mwrfits, obj_strct, objfil, /create
        if (iexp EQ 0) then begin
           allframes = [ obj_strct ]
        endif else begin
           allframes = [ allframes, obj_strct ]
        endelse
     endfor

     ;; ???? At present there is no way to coadd multiple std exposures 
     
     IF NOT KEYWORD_SET(STD) AND NOT KEYWORD_SET(sensfunc) THEN BEGIN
        ;; Combine exposures
        frames = string(objframes.frame)
        
        sxaddpar, master_hdr, 'COMMENT',  'mage_pipe: co-added ' + $
                  strcompress(strjoin(frames, ','), /rem)
        res = 299792.458D/4200.D*objframes[0].slit
        IF KEYWORD_SET(singles) AND N_ELEMENTS(allframes) GT 15 THEN BEGIN
           print, 'Extracting each spectrum before co-add'
           nexp = N_ELEMENTS(allframes)
           FOR kk=0, nexp-1, 15 DO BEGIN
              mage_combspec, allframes[kk:kk+14], fspec
              outflux = 'FSpec/'+strcompress(objframes[0].object,/rem)+ $
                        '_s_'+strtrim(string(kk/15), 2)+'_F.fits'
              outerr  = 'FSpec/'+strcompress(objframes[0].object,/rem)+ $
                        '_s_'+strtrim(string(kk/15), 2)+'_E.fits'
              combname = 'FSpec/'+strcompress(objframes[0].object,/rem)+ $
                         '_s_'+strtrim(string(kk/15), 2)+'_comb.fits'
              mage_1dspec, fspec, outflux, outerr, combname $
                           , hdr=master_hdr, resvel=res 
           ENDFOR

        ENDIF
        ;;Always do the whole combine
        mage_combspec, allframes, fspec
        
        outflux = 'FSpec/'+strcompress(objframes[0].object,/rem)+'_F.fits'
        outerr  = 'FSpec/'+strcompress(objframes[0].object,/rem)+'_E.fits'
        combname = 'FSpec/'+strcompress(objframes[0].object,/rem)+'_comb.fits'
        mage_1dspec, fspec, outflux, outerr, combname $
                     , hdr=master_hdr, resvel=res 
     ENDIF

  ENDFOR

END
