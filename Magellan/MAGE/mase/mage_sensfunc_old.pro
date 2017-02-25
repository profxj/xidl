pro mage_sensfunc, mage, fitsfile, arcfile=arcfile, fluxtable=fluxtable, sensfile=sensfile, pixfile=pixfile, chk=chk, pixflatfile=pixflatfile

  mage_proc, fitsfile, $
             sciimg, scivar, pixflatfile=pixflatfile, hdr=hdr

  ; Generate a 2D skymodel
  print, "Generating the 2D sky model"
  piximg = mrdfits(pixfile)
  orderimg=mrdfits("Orders.fits",0)
  tset_slits = mrdfits("Orders.fits",1)
  skyimage = mage_skymodel(sciimg, scivar, piximg=piximg $
                           , tset_slits=tset_slits)

  ; Generate a wavelenth solution for nearest arc
  print, "Generating arc solution"
  ordr_str=mrdfits("OStr_mage.fits", 1)
  arcfits = arcfile
  mage_arc, arcfits, ordr_str=ordr_str, outarc=arcimg_fil

  ; Find objects.  For bright objects should use the /std flag!
  obj_strct = mage_findobj(sciimg-skyimage, scivar,tset_slits)

  if (keyword_set(CHK)) then begin
     xatv, (sciimg-skyimage) * (orderimg GT 0), min=-20, max=100, /block
;        for qq=0L,14 do xatvplot, obj_strct[qq].trace[0:2048-1], findgen(2048)
  endif

  obj_strct.img_fil = fitsfile
  obj_strct.arc_fil = arcfile

  ; OPTIMALLY EXTRACT
  velpix = 22.0d 
  waveimg = xmrdfits(arcimg_fil)

  m_extechopt, sciimg, sciimg-skyimage, scivar, ordr_str, obj_strct, velpix, /chk, img_arc=waveimg, skyfil=skyimage, helio=0.0, obj_name="junk", ordermask=ordermimage 
  mage_fitstd, fluxtable, obj_strct, sensfile

end 
