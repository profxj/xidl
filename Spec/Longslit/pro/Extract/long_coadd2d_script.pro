
;; Script to create an optimally weighted coadd and then generate a
;; new plan file for input into long_reduce. It then runs
;; long_reduce.pro

scipath = 'Science/'
;; Remove any previous co-adds from the science directory 
spawn, 'rm -f ' + scipath + '*coadd*.fits'
scifiles = findfile(scipath + 'sci-*.fits.gz')
wavefile = findfile('wave*.fits*', count = nfil)
IF KEYWORD_SET(WAVEFILE) THEN waveimg = xmrdfits(wavefile[0], 0) 
IF n_elements(scifiles) EQ 0 THEN message, 'Could not find science files'

;; Output co-add file
outfile = fileandpath(scifiles[0])
outfile = repstr(outfile, 'sci-', 'coadd_2d-')
outfile = repstr(outfile, '.fits.gz', '.fits')

basename = repstr(repstr(fileandpath(scifiles[0]), 'sci-', ''), '.gz', '')
outfile = 'coadd_2d-' + basename 

;; Base S/N weights on the nbright=2 brightest objects
nbright = 2L
weights = long_coadd2d_opt_weights(scifiles, nbright)


long_coadd2d, scifiles, outfile = outfile, weights = weights $
              , imgfinal = imgfinal, subfinal = subfinal $
              , ivarfinal = ivarfinal, maskfinal = maskfinal 
xatv, subfinal, wv = waveimg, /block


basename = 'FORS2_MXU251.1.CHIP1.fits'
outfile = 'coadd_2d-' + basename 
;; Write out a new plan file, with this 
planfile = 'plan_CHIP1.par'
planstr = yanny_readone(planfile, /anonymous)
newplanfile = 'plan_coadd_2d.par'

icals = WHERE(planstr.FLAVOR NE 'science')
icoadd = WHERE(strmatch(planstr.FILENAME, basename))
;;planstr[icoadd].filename = planstr[icoadd].filename ;;outfile
ikeep = [icals, icoadd]
plan_new0 = planstr[ikeep]
nnew = n_elements(plan_new0)
plan_proto = create_struct(name = 'lexp', $
                           'FILENAME', '', $
                           'FLAVOR', '', $
                           'TARGET', '', $
                           'EXPTIME', 0., $
                           'INSTRUMENT', '', $
                           'GRATING', '', $
                           'WAVE', '', $                       
                           'MASKNAME', ''  )
plan_new = replicate(plan_proto, nnew)
;;FOR ii = 0L, nnew-1L DO plan_new[ii].FILENAME =
;;plan_new0[ii].FILENAME
struct_assign, plan_new0, plan_new
yanny_write, newplanfile, ptr_new(plan_new), hdr = planhdr, /align
planstr2 = yanny_readone(newplanfile, hdr = plan, /anonymous)




;;nnew = n_elements(plan_new0)
;;plan_new = replicate(zero_struct(plan_new0[0]), nnew)
;;struct_assign, plan_new0, plan_new

stop
struct_assign, plan_new



FOR ii = 0L, nnew-1L DO struct_assign, plan_new0, plan_new


   stop


indir = yanny_par(planhdr, 'indir', indx = indx)
planhdr[indx] = "indir './' # Raw data diretory"
planhdr = strcompress(planhdr)
yanny_write, newplanfile, ptr_new(plan_new), stnames = 'LEXP', hdr = planhdr, /align

long_reduce, newplanfile

END



