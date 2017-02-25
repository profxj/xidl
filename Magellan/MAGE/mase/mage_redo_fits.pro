arcsave='/Users/jhennawi/mage/mar09_ThAr/Arc1075_old.idl'
restore,arcsave
slit =0.70D
pkwdth=(slit/0.3)*1.2d ;; pixels per slit times 1.2
TOLER = FLOOR(PKWDTH/3.)
mxoff=toler
pksig_vec=[5.0,5.0,replicate(10.0,11),7.0,5.0]
iord=14

pksig=pksig_vec[iord]
path='/Users/jhennawi/mage/mar09_ThAr/calib_arcs_old/'
calibfile=path + 'calib_' + strcompress(string(iord),/rem) + '.sav'
restore,calibfile
incalib=calib_out
calib_out=0
shft=0
stop
x_identify,sv_aspec[*,iord], calib_out, xsize = 1200, ysize = 500 $
           , mfitstr=incalib,mspec=sv_aspec[*,iord] $
           , mshift=shft $
           , linelist =getenv('XIDL_DIR') + '/Spec/Arcs/Lists/mm_thar_mod.lst' $
           , PKSIG =pksig, PKWDTH = pkwdth, TOLER =TOLER, MXOFF = mxoff $
           , /THIN ;;,/FWEIGHT 
stop
newpath='/Users/jhennawi/mage/mar09_ThAr/calib_arcs/'
savefile=newpath + 'calib_' + strcompress(string(iord),/rem) + '.sav'
save,calib_out,file=savefile


;; Why doesn't this work??
;; x_identify,sv_aspec[*,iord], calib_out, xsize = 1200, ysize = 600 $
;;            , linelist =getenv('XIDL_DIR') + '/Spec/Arcs/Lists/mm_thar_mod.lst' $
;;            , PKSIG =psig, PKWDTH = pkwdth, TOLER =TOLER, MXOFF = mxoff $
;;            , /THIN, /FWEIGHT, incalib=all_arcfit[iord]
;stop


END
