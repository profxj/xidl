arcsave='/Users/jhennawi/mage/mar09_ThAr/Arc1075_old.idl'
restore,arcsave
path='/Users/jhennawi/mage/mar09_ThAr/calib_arcs/'
;guessarc = getenv('MAGE_DIR')+'/Calib/MagE_wvguess2.idl'
;restore,guessarc


norders=15L
FOR iord=0L,norders-1L DO BEGIN
   calib_out=0
   calibfile=path + 'calib_' + strcompress(string(iord),/rem) + '.sav'
   restore,calibfile
   all_arcfit[iord]=calib_out
ENDFOR
newguess = getenv('MAGE_DIR')+'/Calib/MagE_wvguess_jfh.idl'
save,all_arcfit,guess_ordr,rejstr,sv_aspec,sv_lines,file=newguess


END
