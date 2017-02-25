
restore,'testbed.sav'
waveimg=mrdfits('Arcs/ArcImg0089.fits.gz')
;;objstr = mage_findobj(sciimg, skyimage, scivar, hdr=hdr,
;;orderfil="OStr_mage.fits",
;;qa='QA/QAFindObj_'+objframes[iexp].fitsfile+'.ps',/CHK)
;;filstd='Object/ObjStr0109.fits'
objstr=mage_findobj(sciimg-skyimage,scivar,tset_slits,filstd=filstd,chk=1)
sciivar=scivar
scihdr=hdr
save,sciimg,sciivar,scihdr,skyimage,piximg,waveimg,tset_slits,objstr,file='ech.sav'

END
