allfiles='Object/' + ['ObjStr1073.fits','ObjStr1074.fits', 'ObjStr2028.fits' $
           , 'ObjStr2030.fits']
;;
calsavfile='GD108_std.sav'
FOR kk=0L,n_elements(allfiles)-1L DO BEGIN
   obj_strct=xmrdfits(allfiles[kk],1)
   mage_flux,calsavfile,obj_strct
   if (kk EQ 0) then allframes = [ obj_strct ] $
   else allframes = [ allframes, obj_strct ]
ENDFOR
mage_combspec, allframes, fspec ;,/chk
mage=mrdfits('mage.fits',1)
imatch=where(strmatch(mage.object,'*SDSS1420*') AND $
             strmatch(mage.exptype,'*SCIENCE*'),nmatch)
master_hdr=headfits(mage[0].RAWPATH + mage[imatch[0]].FITSFILE)
outflux = 'FSpec/SDSS1420bg_F.fits'
outerr  = 'FSpec/SDSS1420bg_E.fits'
combname = 'FSpec/SDSS1420bg_comb.fits'
res = 299792.458/4100.*0.70d
mage_1dspec, fspec, outflux, outerr, combname, hdr=master_hdr, resvel=res

END
