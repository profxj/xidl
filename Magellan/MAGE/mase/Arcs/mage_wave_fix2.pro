
ii=14
igss=14
;guess_fit.LSIG=2.0
;guess_fit.HSIG=2.0

guessarc = getenv('XIDL_DIR')+'/Spec/Longslit/pro/MAGE/ARCHIVE/MagE_wvguess2.idl'
restore,guessarc
;guess_fit=all_arcfit
; x_identify, sv_aspec[*,ii], calib_temp $
;             , mfitstr = guess_fit[igss], mspec =guess_spec[*,igss] $
;             , mshift = shft $
;             , xsize = 1200, ysize = 600 $
;             , linelist = getenv('XIDL_DIR') + $
;             '/Spec/Arcs/Lists/mm_thar.lst' $
;             , PKSIG = psig, PKWDTH = pkwdth, TOLER =TOLER, MXOFF = mxoff $
;             , /THIN, /FWEIGHT 
;stop
;restore,'order6.sav'
;calib_old=calib_temp
;delvarx,calib_temp
psig=5.0
pkwdth=4.0
toler=1.0
mxoff=1.0
x_identify,sv_aspec[*,ii], calib_out, xsize = 1200, ysize = 600 $
            , linelist = getenv('XIDL_DIR') + $
            '/Spec/Arcs/Lists/mm_thar.lst' $
           , PKSIG =psig, PKWDTH = pkwdth, TOLER =TOLER, MXOFF = mxoff $
           , /THIN, /FWEIGHT,incalib=all_arcfit[igss]
stop
x_arclist, linlist, lines
lines.flg_plt = 0
x_templarc, sv_aspec[*,ii], lines, calib_out, MSK=msk $
            , ALL_PK=all_pk, PKWDTH=pkwdth, FORDR=9 $
            ,PKSIG=psig, FLG=flg_templ, BEST=best, /FWEIGHT, $
            /THIN,TOLER=TOLER,mxoff=mxoff
 
gdfit = where(lines.flg_plt EQ 1, ngd)

tmp_fit = {fitstrct}
copy_struct, calib_out, tmp_fit, $
             EXCEPT_TAGS=['FFIT','NRM','RMS']
tmp_fit.flg_rej = 1 
tmp_fit.niter = 3 
tmp_fit.maxrej = 10 
tmp_fit.minpt = 5
          
fit = x_fitrej(lines[gdfit].pix, alog10(lines[gdfit].wave), $
               FITSTR=tmp_fit, REJPT=rejpt)





END
