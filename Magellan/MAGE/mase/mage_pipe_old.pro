pro mage_pipe, mage, obj_id=obj_id, sensfile=sensfile, chk=chk, pixfile=pixfile,  pixflatfile=pixflatfile

  ;setup directories
  IF NOT FILE_TEST('Object', /DIR) THEN FILE_MKDIR, 'Object'
  IF NOT FILE_TEST('Sky', /DIR) THEN FILE_MKDIR, 'Sky'

  ;pick out science targets
  sciframes = mage[where(strtrim(mage.exptype) EQ "SCIENCE" OR strtrim(mage.exptype) EQ "STD"  OR  strtrim(mage.exptype) EQ "BRIGHT" AND mage.obj_id GE 0, nscience)]
  obj_ids = sciframes[uniq(sciframes.obj_id)].obj_id
  
  ;if we only are reducing one object
  if (keyword_set(OBJ_ID)) then begin
     obj_ids = [obj_id]
  endif

  ; Loop over objects
  for iobj=0, n_elements(obj_ids)-1 do begin

     objframes = sciframes[where(sciframes.obj_id EQ obj_ids[iobj], nframes)]
     print, "MAGE_PIPE: Processing object ", objframes[0].object
     
     ;add tags
     master_hdr = objframes[0].hdr
     sxdelpar, master_hdr, "NAXIS1"
     sxdelpar, master_hdr, "CDELT1"
     sxdelpar, master_hdr, "CD1_1"
     sxdelpar, master_hdr, "CRVAL1"
     sxaddpar, master_hdr, "COMMENT", "mage_pipe: v0.2, Build Jan 2009"
     sxaddpar, master_hdr, "COMMENT", "object reduced as: " + objframes[0].exptype
     pipe_time = "mage_pipe: Reduced "+systime()
     sxaddpar, master_hdr, "COMMENT", pipe_time
     
     for iexp=0, nframes-1 do begin
        
        ; Bias subtract and flat field
        mage_proc, objframes[iexp].rawpath+objframes[iexp].fitsfile, $
                   sciimg, scivar, pixflatfile=pixflatfile, hdr=hdr

         ; Generate a 2D skymodel
        print, "Generating the 2D sky model"
        piximg = mrdfits(pixfile)
        orderimg=mrdfits("Orders.fits",0)
        tset_slits = mrdfits("Orders.fits",1)
        skyimage = mage_skymodel(sciimg, scivar, piximg=piximg, $
                                 orderimg=orderimg, tset_slits=tset_slits)
      
        ;Output 2D skymodel
        skyname = 'Sky/'+strtrim(objframes[0].object,2)+'_prelim_sky.fits'
        mwrfits, skyimage, skyname, /create     


        ; Take a look at how we did
        if KEYWORD_SET(CHK) then begin
           print, 'here'
           xatv, (sciimg-skyimage) * (orderimg GT 0), min=-20, max=100, /block
;        for qq=0L,14 do xatvplot, obj_strct[qq].trace[0:2048-1], findgen(2048)
        endif


        ; Generate a wavelenth solution for nearest arc
        print, "Generating arc solution"
        ordr_str=mrdfits("OStr_mage.fits", 1)
        arcfits = objframes[iexp].rawpath+objframes[iexp].arcfile
        ;print, "Arc file: ", arcfits
        mage_arc, arcfits, ordr_str=ordr_str, outarc=arcimg_fil, /clobber

        ; Find objects.  For bright objects should use the /std flag!
        if (objframes[iexp].exptype EQ 'STD' OR objframes[iexp].exptype EQ 'BRIGHT' ) then begin
           obj_strct = mage_findobj(sciimg, skyimage, scivar, hdr=hdr, orderfil="OStr_mage.fits", qa='QA/QAFindObj_'+objframes[iexp].fitsfile+'.ps', /std)
        endif else begin
           obj_strct = mage_findobj(sciimg, skyimage, scivar, hdr=hdr, orderfil="OStr_mage.fits",  qa='QA/QAFindObj_'+objframes[iexp].fitsfile+'.ps')
        endelse

        obj_strct.img_fil = objframes[iexp].rawpath+objframes[iexp].fitsfile
        obj_strct.arc_fil = objframes[iexp].rawpath+objframes[iexp].arcfile

        ; OPTIMALLY EXTRACT
        velpix = 22.0d          ;pixel size in velocity (we should be able to measure this on the image, maybe during the arc process)
        waveimg = xmrdfits(arcimg_fil)
        m_extechopt, sciimg, sciimg-skyimage, scivar, ordr_str, obj_strct, velpix, img_arc=waveimg, skyfil=skyimage, helio=0.0, obj_name=strtrim(objframes[0].object,2), ordermask=ordermimage, model_sky=model_sky, /chk
        ;STOP
        ;Output the final sky model
        skyname = 'Sky/'+strtrim(objframes[0].object,2)+'_final_sky.fits'
        mwrfits, model_sky, skyname, /create    


        mage_flux, sensfile, obj_strct, rej=0.05
        mage_vachelio, obj_strct, head=master_hdr
        
        tmp = strsplit(objframes[iexp].fitsfile, 'mage', /extract)
        outfil = 'Object/ObjStr'+strtrim(tmp[0])
        mwrfits, obj_strct, outfil, /create
        if (iexp EQ 0) then begin
           allframes = [ obj_strct ]
        endif else begin
           allframes = [ allframes, obj_strct ]
        endelse

     endfor

     ; Combine exposures here!

     files = objframes.fitsfile
     msg = "mage_pipe: co-added "+strjoin(files, ', ')
     sxaddpar, master_hdr, "COMMENT", msg
     res = 299792.458/4100. * objframes[0].slit

     mage_combspec, allframes, fspec
s     
     IF NOT FILE_TEST('Final', /DIR) THEN FILE_MKDIR, 'Final'
     
     outflux = 'Final/'+strtrim(objframes[0].object,2)+'_F.fits'
     outerr  = 'Final/'+strtrim(objframes[0].object,2)+'_E.fits'
     combname = 'Final/'+strtrim(objframes[0].object,2)+'_comb.fits'
     
     mage_1dspec, fspec, outflux, outerr, combname, hdr=master_hdr, resvel=res


  endfor

end 
