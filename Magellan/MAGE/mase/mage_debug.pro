allfiles='Object/' + ['ObjStr1062.fits','ObjStr1064.fits' ]
;;
calsavfile='std.sav'
FOR kk=0L,n_elements(allfiles)-1L DO BEGIN
   obj_strct=xmrdfits(allfiles[kk],1)
   mage_flux,calsavfile,obj_strct ;,/chk
   if (kk EQ 0) then allframes = [ obj_strct ] $
   else allframes = [ allframes, obj_strct ]
ENDFOR
mage_combspec, allframes, fspec,/chk
mage=mrdfits('mage_1144.fits',1)
imatch=where(strmatch(mage.object,'*SDSS1144*') AND $
             strmatch(mage.exptype,'*SCIENCE*'),nmatch)
master_hdr=headfits(mage[0].RAWPATH + mage[imatch[0]].FITSFILE)
outflux = 'FSpec/SDSS1144fg_test_F.fits'
outerr  = 'FSpec/SDSS1144fg_test_E.fits'
combname = 'FSpec/SDSS1144fg_test_comb.fits'
res = 299792.458/4100.*0.70d
mage_1dspec, fspec, outflux, outerr, combname, hdr=master_hdr, resvel=res

END
