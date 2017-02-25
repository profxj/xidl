;; maskpath = path to masks, i.e. ~/redux/Mask_Name/Blue600/'
;; iref = reference image, that long_coadd will scale things to for
;; determining SNR
;; CHECK = interactive
;; NBRIGHT = number of bright objects to average over for weight determination
PRO LONG_COADD2D_RUN, maskpath, iref = iref, CHECK = CHECK, NBRIGHT = NBRIGHT1

IF n_elements(iref) EQ 0 THEN iref = 0
IF KEYWORD_SET(NBRIGHT1) THEN NBRIGHT = NBRIGHT1 ELSE NBRIGHT = 1
;; Script to create an optimally weighted coadd and then generate a
;; new plan file for input into long_reduce.   
;; Remove any previous co-adds from the science directory 
spawn, 'rm -f ' + maskpath + '*coadd*.fits'
scifiles = findfile(maskpath + 'Science/sci-*.fits.gz')
wavefile = findfile(maskpath + 'wave*.fits*', count = nfil)
IF KEYWORD_SET(WAVEFILE) THEN waveimg = xmrdfits(wavefile[0], 0) 
IF n_elements(scifiles) EQ 0 THEN message, 'Could not find science files'

;; Output co-add file
outfile = fileandpath(scifiles[0])
outfile = repstr(outfile, 'sci-', 'coadd_2d-')
outfile = maskpath + repstr(outfile, '.fits.gz', '.fits')

basename = repstr(repstr(fileandpath(scifiles[0]), 'sci-', ''), '.gz', '')
outfile = maskpath + 'coadd_2d-' + basename 
;; Base S/N weights on the nbright=1 brightest objects
weights = long_coadd2d_opt_weights(scifiles, nbright, iref = iref $
                                   , ymult = ymult)
long_coadd2d, scifiles, outfile = outfile, weights = weights $
              , imgfinal = imgfinal, subfinal = subfinal $
              , ivarfinal = ivarfinal, maskfinal = maskfinal  $
              , ymult = ymult, waveimg = waveimg

IF KEYWORD_SET(CHECK) THEN BEGIN
   xatv, subfinal*maskfinal, wv = waveimg ;;, /block
   objstruct = xmrdfits(scifiles[iref], 5)
   xatvplot, objstruct.xpos, objstruct.ypos, psym = 3
ENDIF

END
