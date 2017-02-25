;+
; NAME:
;   tspec_diff_proc
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;-
PRO TSPEC_DIFF_PROC, filenames, pixflatfile, illumflatfile $
                     , tset_slits, ZAP = ZAP $
                     , AB = AB, ivar_AB = ivar_AB, piximg = piximg $
                     , targdir = targdir, sky_AB = sky_AB $
                     , waveimg = waveimg, TELLURIC = TELLURIC, hdr = hdr_A $
                     , WVCHK = WVCHK
;, ivar_pos = ivar_pos, ivar_neg = ivar_neg $
IF NOT KEYWORD_SET(scatt_shift) THEN scatt_shift = 5L
;; This cuts will act like a bad pixel mask
minval = 0.7
maxval = 1.3
IF NOT KEYWORD_SET(BSP) THEN BSP = 0.5
IF NOT KEYWORD_SET(CR_THRESH) THEN CR_THRESH = 10.0
IF NOT KEYWORD_SET(CR_VAL) THEN CR_VAL    = 100.0

dimt = size(tset_slits.coeff, /dimen)
ordermask = tspec_ordermask(tset_slits, order_vec = order_vec)
norders = dimt[1]
ximg = long_slits2x(tset_slits, edgmask = edgmask)
dims = size(ordermask, /dimen)
nx = dims[0]
ny = dims[1]

; Read in science object files 
tspec_proc, filenames[0], img_A, ivar_A, hdr = hdr_A $
            , pixflatfile = pixflatfile, illumflatfile = illumflatfile
tspec_proc, filenames[1], img_B, ivar_B, hdr = hdr_B $
            , pixflatfile = pixflatfile, illumflatfile = illumflatfile
maskstack = lonarr(nx*ny, 2)
img_stack = fltarr(nx*ny, 2)
img_stack[*, 0] = img_A
img_stack[*, 1] = img_B
avg_sky =  djs_avsigclip(img_stack, 2, sigrej = 2.0, inmask = maskstack $
                         , outmask = outmask)
avg_sky = reform(avg_sky, nx, ny)
piximg = tspec_makescipix(avg_sky, tset_slits, pixset = pixset, chk = wvchk)

;; KLUDGE FOR NOW
;waveimg = tspec_waveimg(piximg, ordermask, order_vec)

i0 = strpos(filenames[0],'/', /reverse_sear)
i1 = strpos(filenames[0],'.fits', /reverse_sear) - 1
qafile = 'QA/'+strmid(filenames[0],i0+1,i1-i0)+'_qawave.ps'
waveimg = tspec_waveimg(avg_sky, hdr_A, tset_slits, piximg = piximg $
                        , QAFILE = QAFILE, CHK = WVCHK)

;gnirs_reject_cr, a1, ivar_a1, a2, ivar_a2, tset_slits, piximg, ordermask $
;                 , cr_mask1 = cr_mask_a1, cr_mask2 = cr_mask_a2
;gnirs_reject_cr, b1, ivar_b1, b2, ivar_b2, tset_slits, piximg, ordermask $
;                 , cr_mask1 = cr_mask_b1, cr_mask2 = cr_mask_b2
;maskstack[*, 0] = (cr_mask_a1 EQ 0)
;maskstack[*, 1] = (cr_mask_a2 EQ 0)
;maskstack[*, 2] = (cr_mask_b1 EQ 0)
;maskstack[*, 3] = (cr_mask_b2 EQ 0)
;avg_sky =  djs_avsigclip(img_stack, 2, sigrej = 1.3, inmask = maskstack $
;                         , outmask = outmask)
;wpix_img = reform(avg_sky, nx, ny)
;print, "Wavelength calibration"
; Generate final wavelength calibration wavelength map 
;waveimg = gnirs_waveimg(wpix_img, wpixhdr, tset_slits, piximg = piximg $
;                        , QAFILE = QAFILE, CHK = WVCHK)
IF KEYWORD_SET(TELLURIC) THEN BEGIN
;   Read in Telluric files now
; Read in science object files 
tspec_proc, telluric[0], img_A, ivar_A, hdr = hdr_A $
            , pixflatfile = pixflatfile, illumflatfile = illumflatfile
tspec_proc, telluric[1], img_B, ivar_B, hdr = hdr_B $
            , pixflatfile = pixflatfile, illumflatfile = illumflatfile
maskstack = lonarr(nx*ny, 2)
img_stack = fltarr(nx*ny, 2)
img_stack[*, 0] = img_A
img_stack[*, 1] = img_B
avg_sky =  djs_avsigclip(img_stack, 2, sigrej = 2.0, inmask = maskstack $
                         , outmask = outmask)
avg_sky = reform(avg_sky, nx, ny)
ENDIF

;;   First pass AB sky subtraction.  Also corrects for bias and dark counts
print, "  First pass AB sky subtraction"
AB = img_A-img_B
sig2_AB = (ivar_A GT 0)/(ivar_A + (ivar_A LE 0.0)) + $
          (ivar_B GT 0)/(ivar_B + (ivar_B LE 0.0)) 
ivar_AB = (sig2_AB GT 0.0)/(sig2_AB + (sig2_AB LE 0))
finalmask = (ivar_A GT 0) AND (ivar_B GT 0)
sky_AB  = avg_sky

chk = 1
;;   Bspline sky subtraction for CR rejection
IF KEYWORD_SET(TELLURIC) THEN BEGIN
;   Skip second-pass for standard stars and bright objects.
   print, "    Bright calibration source; skipping 2nd pass..."
   sky_model = sky_AB 
   sky_resids = 0.0*sky_AB
ENDIF ELSE BEGIN
;      Otherwise do 2nd pass
   print, "  CR sky subtraction..."
   y_img = replicate(1.0, nx)#findgen(ny)
   sky_resids = img_A*0.0
   FOR iorder = 0, norders-1 DO BEGIN
      inorder = where(ordermask EQ order_vec[iorder], nord)
      fitpix = where(ordermask EQ order_vec[iorder] AND $
                     finite(img_A)  AND $
                     edgmask EQ 0 AND $
                     abs(img_A) LE 1.0d6 AND $
                     ivar_A GT 0.0, npix)
;                       ximg GT 0.1  AND $
;                       ximg LT 0.9  AND $
      psort = sort(piximg[fitpix])
      sset = bspline_iterfit(piximg[fitpix[psort]], img_A[fitpix[psort]] $
                             , invvar = (ivar_A[fitpix[psort]] GT 0.0) $
                             , upper = 3, lower = 3 $
                             , bkspace = bsp, maxiter = 20, maxrej = 10 $
                             , /silent)
      sky_resids[inorder] = bspline_valu(piximg[inorder], sset)
      sky_model = sky_AB + sky_resids
      IF KEYWORD_SET(CHK) THEN BEGIN
         plotx = piximg[fitpix[psort]]
         ploty = img_A[fitpix[psort]]
         rms = sqrt(djs_median(ploty^2))
         ;x_splot, plotx, ploty, psym1 = 3 $
         ;        , ymnx = [-10.0*rms, 10.0*rms] $
         ;        , xtwo = plotx, ytwo = bspline_valu(plotx, sset), /block
      ENDIF
   ENDFOR
ENDELSE

;; Reject CRs, first pass
FWHM = 4.0                  
sigma_psf = FWHM/2.35482D
psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
cr_sharp_A = psf_reject_cr(AB-sky_resids, ivar_AB, psfvals $
                           , satmask = ((AB-sky_resids) LT -10))
cr_sharp_B = psf_reject_cr(sky_resids-AB, ivar_AB, psfvals $
                           , satmask = ((sky_resids-AB) LT -10))

ivar_A = (cr_sharp_A EQ 0)*ivar_A
ivar_B = (cr_sharp_B EQ 0)*ivar_B
;;   First pass AB sky subtraction.  Also corrects for bias and dark counts
print, "  First pass AB sky subtraction"
AB = img_A-img_B
sig2_AB = (ivar_A GT 0)/(ivar_A + (ivar_A LE 0.0)) + $
          (ivar_B GT 0)/(ivar_B + (ivar_B LE 0.0)) 
ivar_AB = (sig2_AB GT 0.0)/(sig2_AB + (sig2_AB LE 0))
finalmask = (ivar_A GT 0) 
sky_AB  = avg_sky


cr_mask = (cr_sharp_A EQ 1) OR (cr_sharp_B EQ 1)
cr_mask = long(cr_mask EQ 0)
;; Reject CRs, first pass
;FWHM = 3.0                  
;sigma_psf = FWHM/2.35482D
;psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
;cr_sharp_A = psf_reject_cr(AB, ivar_AB, psfvals $
;                           , satmask = (AB LT -10))
;cr_sharp_B = psf_reject_cr(-AB, ivar_AB, psfvals $
;                           , satmask = (-AB LT -10))
;X2_AB = (AB-sky_resids)*sqrt(ivar_AB)

;   Mask out-of-order pixels
out_order = where(ordermask EQ 0)
ivar_AB[out_order] = 0.0
finalmask[out_order] = 0.0

RETURN
END
