PRO GNIRS_ABBA_PROC, filenames1, flat, tset_slits, acqfiles, ZAP = ZAP $
                     , abba = abba, ivar_abba = ivar_abba, piximg = piximg $
                     , targdir = targdir, sky_abba = sky_abba $
                     , waveimg = waveimg, TELLURIC = TELLURIC, hdr = hdr1 $
                     , WVCHK = WVCHK
;, ivar_pos = ivar_pos, ivar_neg = ivar_neg $
IF NOT KEYWORD_SET(scatt_shift) THEN scatt_shift = 5L
;; This cuts will act like a bad pixel mask
minval = 0.7
maxval = 1.3

dims = size(flat, /dimen)
nx = dims[0]
ny = dims[1]

ordermask = long_slits2mask(tset_slits)
ordermask[WHERE(ordermask GT 0)] = ordermask[WHERE(ordermask GT 0)] + 2L
ximg = long_slits2x(tset_slits)

IF KEYWORD_SET(TELLURIC) THEN filenames = telluric $
ELSE filenames = filenames1

; Read in science object files 
gnirs_proc, filenames[0], a1, ivar_a1, hdr = hdr1, flatimg = flat
gnirs_proc, filenames[1], b1, ivar_b1, hdr = hdr2, flatimg = flat
gnirs_proc, filenames[2], b2, ivar_b2, hdr = hdr3, flatimg = flat
gnirs_proc, filenames[3], a2, ivar_a2, hdr = hdr4, flatimg = flat
;    Setting equal to the third in an ABBA sequence sets the UT start of
;    the exposure to somewhere near the midpoint of the sequence.
;    However, we do want to keep the exposure time equal to a single
;    exposure rather than adding the 4.  This is because the exposures
;    are normalized to counts = 0.5 * (A1 + A2 - B1 - B2)
hdr = hdr3
badpixmask = (flat GT minval AND flat LT maxval)
badpixmask[*, 1012:*] = 0      ; zero out y > 1012 bad regions on chip
ivar_a1 = ivar_a1*badpixmask
ivar_a2 = ivar_a2*badpixmask
ivar_b1 = ivar_b1*badpixmask
ivar_b2 = ivar_b2*badpixmask
;; No CR zapping, pattern noise removal, or scattered light/persistence
;; removal for Tellurics

;; Remove pattern noise
splog, "  Removing GNIRS pattern noise"
gnirs_remove_patternnoise, a1
gnirs_remove_patternnoise, a2
gnirs_remove_patternnoise, b1
gnirs_remove_patternnoise, b2

;gnirs_reject_cr, a1, ivar_a1, a2, ivar_a2, tset_slits, piximg, ordermask $
;                 , cr_mask1 = cr_mask_a1, cr_mask2 = cr_mask_a2
;gnirs_reject_cr, b1, ivar_b1, b2, ivar_b2, tset_slits, ordermask $
;                 , cr_mask1 = cr_mask_b1, cr_mask2 = cr_mask_b2
;; Remove scattered light and persistence, read in acquisition files
nacq = n_elements(acqfiles)
acq_stk = fltarr(nx, ny, nacq)
FOR j = 0L, nacq -1L DO BEGIN
   gnirs_proc, acqfiles[j], acq1, hdr = hdr, /pattern
   acq_stk[*, *, j] = acq1
ENDFOR
;; Creat template for acquisition persistence pattern
acq_bar = djs_avsigclip(acq_stk, 3) 
djs_iterstat, acq_bar[WHERE(acq_bar GT 1000.0)], mean = mean, median = median
acq_template = double(acq_bar GT mean/5.0D)
;; Create a model image of the acquisition field persistence pattern.
acq_template[*, 510:*] = 0
acq_template[*, 0:330] = 0
acq_template[0:160, *] = 0
acq_template[880:*, *] = 0
;; use reject_cr to fill in sharp features in the acq_template
FWHM = 3.0                  
sigma_psf = FWHM/2.35482D
psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
acq_temp_ivar = fltarr(nx, ny) + 1.0/(0.01)^2
acq_cr = psf_reject_cr((NOT acq_template), acq_temp_ivar, psfvals $
                       , satmask = (acq_template GT 8.0d4))
acq_template = acq_template OR acq_cr
;; Correct for zero point amplifier offsets between quadrants
xval = findgen(nx) # replicate(1.0, ny)
yval = replicate(1.0, nx) # findgen(ny)
;;  quadrants
ul = xval GE 0    AND xval LE nx/2 AND yval GE ny/2 AND yval LE ny-1L
ur = xval GE nx/2 AND xval LE nx-1 AND yval GE ny/2 AND yval LE ny-1L
ll = xval GE 0    AND xval LE nx/2 AND yval GE 0    AND yval LE ny/2
lr = xval GE nx/2 AND xval LE nx-1 AND yval GE 0    AND yval LE ny/2
;; remove offset from a1
ul_ind = WHERE(ul AND ordermask EQ 0 AND acq_template EQ 0)
ur_ind = WHERE(ur AND ordermask EQ 0 AND acq_template EQ 0)
ll_ind = WHERE(ll AND ordermask EQ 0 AND acq_template EQ 0)
lr_ind = WHERE(lr AND ordermask EQ 0 AND acq_template EQ 0)
djs_iterstat, a1[ul_ind], mean = mean_ul, median = med_ul, sigrej = 2.0
a1[WHERE(ul)] -=  med_ul
djs_iterstat, a1[ur_ind], mean = mean_ur, median = med_ur, sigrej = 2.0
a1[WHERE(ur)] -=  med_ur 
djs_iterstat, a1[ll_ind], mean = mean_ll, median = med_ll, sigrej = 2.0
a1[WHERE(ll)] -=  med_ll 
djs_iterstat, a1[lr_ind], mean = mean_lr, median = med_lr, sigrej = 2.0
a1[WHERE(lr)] -=  med_lr
;; remove offset from a2
ul_ind = WHERE(ul AND ordermask EQ 0 AND acq_template EQ 0)
ur_ind = WHERE(ur AND ordermask EQ 0 AND acq_template EQ 0)
ll_ind = WHERE(ll AND ordermask EQ 0 AND acq_template EQ 0)
lr_ind = WHERE(lr AND ordermask EQ 0 AND acq_template EQ 0)
djs_iterstat, a2[ul_ind], mean = mean_ul, median = med_ul, sigrej = 2.0
a2[WHERE(ul)] -=  med_ul
djs_iterstat, a2[ur_ind], mean = mean_ur, median = med_ur, sigrej = 2.0
a2[WHERE(ur)] -=  med_ur 
djs_iterstat, a2[ll_ind], mean = mean_ll, median = med_ll, sigrej = 2.0
a2[WHERE(ll)] -=  med_ll 
djs_iterstat, a2[lr_ind], mean = mean_lr, median = med_lr, sigrej = 2.0
a2[WHERE(lr)] -=  med_lr
;; remove offset from b1
ul_ind = WHERE(ul AND ordermask EQ 0 AND acq_template EQ 0)
ur_ind = WHERE(ur AND ordermask EQ 0 AND acq_template EQ 0)
ll_ind = WHERE(ll AND ordermask EQ 0 AND acq_template EQ 0)
lr_ind = WHERE(lr AND ordermask EQ 0 AND acq_template EQ 0)
djs_iterstat, b1[ul_ind], mean = mean_ul, median = med_ul, sigrej = 2.0
b1[WHERE(ul)] -=  med_ul
djs_iterstat, b1[ur_ind], mean = mean_ur, median = med_ur, sigrej = 2.0
b1[WHERE(ur)] -=  med_ur 
djs_iterstat, b1[ll_ind], mean = mean_ll, median = med_ll, sigrej = 2.0
b1[WHERE(ll)] -=  med_ll 
djs_iterstat, b1[lr_ind], mean = mean_lr, median = med_lr, sigrej = 2.0
b1[WHERE(lr)] -=  med_lr
;; remove offset from b2
ul_ind = WHERE(ul AND ordermask EQ 0 AND acq_template EQ 0)
ur_ind = WHERE(ur AND ordermask EQ 0 AND acq_template EQ 0)
ll_ind = WHERE(ll AND ordermask EQ 0 AND acq_template EQ 0)
lr_ind = WHERE(lr AND ordermask EQ 0 AND acq_template EQ 0)
djs_iterstat, b2[ul_ind], mean = mean_ul, median = med_ul, sigrej = 2.0
b2[WHERE(ul)] -=  med_ul
djs_iterstat, b2[ur_ind], mean = mean_ur, median = med_ur, sigrej = 2.0
b2[WHERE(ur)] -=  med_ur 
djs_iterstat, b2[ll_ind], mean = mean_ll, median = med_ll, sigrej = 2.0
b2[WHERE(ll)] -=  med_ll 
djs_iterstat, b2[lr_ind], mean = mean_lr, median = med_lr, sigrej = 2.0
b2[WHERE(lr)] -=  med_lr

splog, 'Removing GNIRS persistence/scatt light'

slitmask_left  = long_slits2mask(tset_slits, xshift = -abs(scatt_shift))
slitmask_right = long_slits2mask(tset_slits, xshift = abs(scatt_shift))
scattmask = (slitmask_left EQ 0 AND slitmask_right EQ 0)
notedg = (smooth(acq_template, [2, 2]) - acq_template) EQ 0

itemp  = WHERE(notedg AND acq_template AND scattmask, ntemp)
djs_iterstat, a1[itemp], sigrej = 2.0, mean = mean_a1
a1 = a1 - mean_a1*acq_template
itemp  = WHERE(notedg AND acq_template AND scattmask, ntemp)
djs_iterstat, a2[itemp], sigrej = 2.0, mean = mean_a2
a2 = a2 - mean_a2*acq_template
itemp  = WHERE(notedg AND acq_template AND scattmask, ntemp)
djs_iterstat, b1[itemp], sigrej = 2.0, mean = mean_b1
b1 = b1 - mean_b1*acq_template
itemp  = WHERE(notedg AND acq_template AND scattmask, ntemp)
djs_iterstat, b2[itemp], sigrej = 2.0, mean = mean_b2
b2 = b2 - mean_b2*acq_template

;; First pass wavelength calibration to aid CR rejection
split1 = strsplit(filenames[3], 'S', /extract)
nsplit = n_elements(split1)
split2 = strsplit(split1[nsplit-1L], '.fits*', /extract)
qafile = targdir + '/waveQA-' + gnirs_fileprefix(filenames[0]) $
  + '-' + strcompress(string(long(split2[0])), /REM) + '.ps'
wpixhdr = hdr1
maskstack = lonarr(nx*ny, 4)
img_stack = fltarr(nx*ny, 4)
img_stack[*, 0] = a1
img_stack[*, 1] = a2
img_stack[*, 2] = b1
img_stack[*, 3] = b2
maskstack[*, 0] = (ivar_a1 LE 0)
maskstack[*, 1] = (ivar_a2 LE 0)
maskstack[*, 2] = (ivar_b1 LE 0)
maskstack[*, 3] = (ivar_b2 LE 0)
avg_sky =  djs_avsigclip(img_stack, 2, sigrej = 1.3, inmask = maskstack $
                         , outmask = outmask)
wpix_img = reform(avg_sky, nx, ny)
waveimg = gnirs_waveimg(wpix_img, wpixhdr, tset_slits, piximg = piximg $
                        , QAFILE = QAFILE, CHK = WVCHK)
gnirs_reject_cr, a1, ivar_a1, a2, ivar_a2, tset_slits, piximg, ordermask $
                 , cr_mask1 = cr_mask_a1, cr_mask2 = cr_mask_a2
gnirs_reject_cr, b1, ivar_b1, b2, ivar_b2, tset_slits, piximg, ordermask $
                 , cr_mask1 = cr_mask_b1, cr_mask2 = cr_mask_b2
maskstack[*, 0] = (cr_mask_a1 EQ 0)
maskstack[*, 1] = (cr_mask_a2 EQ 0)
maskstack[*, 2] = (cr_mask_b1 EQ 0)
maskstack[*, 3] = (cr_mask_b2 EQ 0)
avg_sky =  djs_avsigclip(img_stack, 2, sigrej = 1.3, inmask = maskstack $
                         , outmask = outmask)
wpix_img = reform(avg_sky, nx, ny)
print, "Wavelength calibration"
; Generate final wavelength calibration wavelength map 
waveimg = gnirs_waveimg(wpix_img, wpixhdr, tset_slits, piximg = piximg $
                        , QAFILE = QAFILE, CHK = WVCHK)
IF KEYWORD_SET(TELLURIC) THEN BEGIN
;   Read in Telluric files now
    gnirs_proc, filenames1[0], a1, ivar_a1, hdr = hdr1, /pattern
    gnirs_proc, filenames1[1], b1, ivar_b1, hdr = hdr2, /pattern
    gnirs_proc, filenames1[2], b2, ivar_b2, hdr = hdr3, /pattern
    gnirs_proc, filenames1[3], a2, ivar_a2, hdr = hdr4, /pattern
    gnirs_reject_cr, a1, ivar_a1, a2, ivar_a2, tset_slits, piximg, ordermask $
                     , cr_mask1 = cr_mask_a1, cr_mask2 = cr_mask_a2
    gnirs_reject_cr, b1, ivar_b1, b2, ivar_b2, tset_slits, piximg, ordermask $
                     , cr_mask1 = cr_mask_b1, cr_mask2 = cr_mask_b2
ENDIF

; average the A and B images using the CR masks
amask_sum = cr_mask_a1 + cr_mask_a2
amask = cr_mask_a1 OR cr_mask_a2
abar = (cr_mask_a1*a1 + cr_mask_a2*a2)/(amask_sum + (amask_sum EQ 0))
sig2a_bar =  (cr_mask_a1/(ivar_a1 + (ivar_a1 EQ 0)) $
              +  cr_mask_a2/(ivar_a2 + (ivar_a2 EQ 0))) $
  /(amask_sum^2 + (amask_sum EQ 0))
;ivar_pos = amask*amask_sum^2/(sig2a_sum + (sig2a_sum EQ 0))

bmask_sum = cr_mask_b1 + cr_mask_b2
bmask = cr_mask_b1 OR cr_mask_b2
bbar = (cr_mask_b1*b1 + cr_mask_b2*b2)/(bmask_sum + (bmask_sum EQ 0))
sig2b_bar =  (cr_mask_b1/(ivar_b1 + (ivar_b1 EQ 0)) $
              +  cr_mask_b2/(ivar_b2 + (ivar_b2 EQ 0))) $
  /(bmask_sum^2 + (bmask_sum EQ 0))
;ivar_neg = bmask*bmask_sum^2/(sig2b_sum + (sig2b_sum EQ 0))
finalmask = badpixmask*amask*bmask
;   First pass ABBA sky subtraction.  Also corrects for bias and dark counts
print, "  First pass ABBA sky subtraction"
abba = abar-bbar
sig2_abba = sig2a_bar + sig2b_bar
ivar_abba = finalmask*(sig2_abba GT 0.0)/(sig2_abba + (sig2_abba EQ 0))

; ABBA inverse variance is needed for bspline fitting. Since we don't 
; want the objects in the ABBA bspline we use sky only noise for both 
; making the object pixels more likely to be rejected 
;; THIS IS ALL OUTDATED SINCE WE ARE NOT USING THE S/N WEIGHTS IN THE SKY
;; SUBTRACTION BECAUSE OF TOO MANY HOT PIXELS
;ivar_abba = (ordermask GT 0.0 AND ximg LE 0.5)*ivar_neg + $
;  (ordermask GT 0.0 AND ximg GT 0.5)*ivar_pos
sky_abba  = (ordermask GT 0.0 AND ximg LE 0.5)*bbar + $
  (ordermask GT 0.0 AND ximg GT 0.5)*abar
;   Mask out-of-order pixels
out_order = where(ordermask EQ 0)
;abba[out_order] = 0
;ivar_pos[out_order] = 0.0
;ivar_neg[out_order] = 0.0
ivar_abba[out_order] = 0.0
finalmask[out_order] = 0.0
; Mask all with combined mask
;ivar_pos  = ivar_pos*finalmask
;ivar_neg  = ivar_neg*finalmask


;; Edge artifacts are still appearing in the abba image. I think this is 
;; partly due to variations in how the slit is illimuninated as a function 
;; of time position angle but I can't track down exactly why. This will 
;; genereally be a problem for any differencing though I suspect but I don't
;; understand what is making the edges different.???

RETURN
END
