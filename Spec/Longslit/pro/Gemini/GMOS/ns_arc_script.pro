
slitfile = 'slitmask.fits'
;filename = 'N20080407S0192.fits.gz'
;wavefile = 'wave-N20080407S0192-ccd3.fits'
;qafile =  'wave-N20080407S0192-ccd3.ps'
;savefile = repstr(wavefile, '.fits', '.sav')

filename = 'N20080407S0194.fits.gz'
wavefile = 'wave-N20080407S0194-ccd2.fits'
qafile =  'wave-N20080407S0194-ccd2.ps'
savefile = repstr(wavefile, '.fits', '.sav')

slit_arr = mrdfits(slitfile, 0)
tset_slits1 = mrdfits(slitfile, 1)
tset_slits2 = mrdfits(slitfile, 2)
tset_slits3 = mrdfits(slitfile, 3)

;; expand slit edges
; traceset2xy, tset_slits1[0], rows1, left1
; traceset2xy, tset_slits1[1], rows1, right1
; traceset2xy, tset_slits2[0], rows2, left2
; traceset2xy, tset_slits2[1], rows2, right2
; traceset2xy, tset_slits3[0], rows3, left3
; traceset2xy, tset_slits3[1], rows3, right3

; if (size(left1, /n_dimen) EQ 1) then nslit1 = 1 $
; else nslit1 = (size(left1, /dimens))[1]
; if (size(left2, /n_dimen) EQ 1) then nslit2 = 1 $
; else nslit2 = (size(left2, /dimens))[1]
; if (size(left3, /n_dimen) EQ 1) then nslit3 = 1 $
; else nslit3 = (size(left3, /dimens))[1]

;;----------
;; Read the arc image
long_proc, filename, arcimg, arcivar, hdr = hdr, ccdonly = 2

wstruct = long_wstruct(hdr[*, 0], linelist = linelist, reid_file = reid_file)

xfit = long_waveimg(arcimg, arcivar, tset_slits2, wstruct, savefile $
                    , fwhmset = fwhmset, qafile = qafile, box_rad = 2 $
                    , FWCOEFF = 2)

pixset = long_wavepix(arcimg, tset_slits2, fwhm = fwhmset.MEDIAN $
                      , box_radius = wstruct.radius $
                      , sig_thresh = wstruct.sig_wpix $
                      , pkwdth = wstruct.pkwdth $
                      , TOLER = wstruct.TOLER, CHK = CHK)

piximg = long_wpix2image(pixset, tset_slits2, XFIT = xfit $
                         , waveimg = waveimg)


END
