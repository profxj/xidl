arcsave='/Users/jhennawi/mage/mar09_ThAr/Arc1075.idl'
restore,arcsave
iord=10
slit =0.70D
pkwdth=round(slit/1.0d*6L)
TOLER = FLOOR(PKWDTH/3.)
mxoff=toler
psig=20.0

;; Why doesn't this work??
x_identify,sv_aspec[*,iord], calib_out, xsize = 1200, ysize = 600 $
           , linelist = getenv('XIDL_DIR') + '/Spec/Arcs/Lists/mm_thar1.lst' $
           , PKSIG =psig, PKWDTH = pkwdth, TOLER =TOLER, MXOFF = mxoff $
           , /THIN, /FWEIGHT, incalib=all_arcfit[iord]
stop
x_identify,sv_aspec[*,iord], calib_out, xsize = 1200, ysize = 600 $
           ,mfitstr=all_arcfit[iord],mspec=sv_aspec[*,iord] $
           ,mshift=shft $
           , linelist = getenv('XIDL_DIR') + '/Spec/Arcs/Lists/mm_thar1.lst' $
           , PKSIG =psig, PKWDTH = pkwdth, TOLER =TOLER, MXOFF = mxoff $
           , /THIN, /FWEIGHT ;;;, incalib=all_arcfit[iord]



stop
path='/Users/jhennawi/mage/mar09_ThAr/calib_arcs/'
savefile=path + 'calib_' + strcompress(string(iord),/rem) + '.sav'
save,calib_out,file=savefile



END
