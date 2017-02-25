PRO gnirs_reject_cr, a1, ivar_a1, a2, ivar_a2, tset_slits, piximg, ordermask $
                     , cr_mask1 = cr_mask1, cr_mask2 = cr_mask2 $
                     , CR_THRESH = CR_THRESH, CR_VAL = CR_VAL

IF NOT KEYWORD_SET(BSP) THEN BSP = 0.8
IF NOT KEYWORD_SET(CR_THRESH) THEN CR_THRESH = 10.0
IF NOT KEYWORD_SET(CR_VAL) THEN CR_VAL    = 1000.0

ximg = long_slits2x(tset_slits, edgmask = edgmask)

order_vec = [3, 4, 5, 6, 7, 8]
dims = size(a1, /dimens)
nx = dims[0]
ny = dims[1]
dimt = size(tset_slits.coeff, /dimen)
norders = dimt[1]
; generate left and right edge of slits
;traceset2xy, tset_slits[0], rows, left_edge
;traceset2xy, tset_slits[1], rows, right_edge
;trace = (left_edge + right_edge)/2.0D
; Median filtering is more robust against cosmics
;spec_a1 = fltarr(ny, norders)
;spec_a2 = fltarr(ny, norders)
;BOX_RAD = 10
;FOR iorder = 0L, norders-1L DO BEGIN
;    FOR j = 0L, ny-1L DO BEGIN
;        left  = floor(trace[j, iorder] - BOX_RAD)
;        right = ceil(trace[j, iorder] + BOX_RAD)
;        sub_a1 = a1[left:right, j]
;        sub_a2 = a2[left:right, j]
;        djs_iterstat, sub_a1, median = median1, sigrej = 2.0
;        spec_a1[j, iorder] = median1
;        djs_iterstat, sub_a2, median = median2, sigrej = 2.0
;        spec_a2[j, iorder] = median2
;    ENDFOR
;ENDFOR

;npoly = 4
;nback = 1
;scale_12 = fltarr(ny, norders)
;yvector = findgen(ny)/double(ny-1)
sky_resids = fltarr(nx, ny)
;addimg = fltarr(nx, ny)
y_img = findgen(ny)## replicate(1.0, nx)

adiff_cr = (a1-a2)
var_a1   = (ivar_a1 GT 0.0)/(ivar_a1 + (ivar_a1 LE 0.0))
var_a2   = (ivar_a2 GT 0.0)/(ivar_a2 + (ivar_a2 LE 0.0))
var_adiff  = var_a1 + var_a2
ivar_adiff = (var_adiff GT 0.0)/(var_adiff + (var_adiff LE 0.0))

FOR iorder = 0L, norders-1L DO BEGIN
    inorder = where(ordermask EQ order_vec[iorder], nord)
    fitpix = where(ordermask EQ order_vec[iorder] AND $
                   finite(adiff_cr)  AND $
                   edgmask EQ 0 AND $
                   abs(adiff_cr) LE 5.0d4 AND $
                   ivar_adiff GT 0.0, npix)
    psort = sort(piximg[fitpix])
    sset = bspline_iterfit(piximg[fitpix[psort]], adiff_cr[fitpix[psort]] $
                           , invvar = (ivar_adiff[fitpix[psort]] GT 0.0) $
                           , upper = 3, lower = 3 $
                           , bkspace = bsp, maxiter = 20, maxrej = 10 $
                           , /silent)
    sky_resids[inorder] = bspline_valu(piximg[inorder], sset)
ENDFOR
adiff_final = adiff_cr - sky_resids
;   Take the PSF width to be that of the spectra direction (3 pixels) for 
;   the 3-pixel slit. This prevents the routine from rejecting sky lines
;   This is a description of the 3x3 core of the 2D PSF for reject_cr.pro
;    
;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
;    PSFVALS[0]          1.   PSFVALS[0]
;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
FWHM = 3.0                  
sigma_psf = FWHM/2.35482D
psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
cr_sharp1 = psf_reject_cr(adiff_final, ivar_a1, psfvals $
                           , satmask = (adiff_final LT -10))
cr_sharp2 = psf_reject_cr(-adiff_final, ivar_a2, psfvals $
                          , satmask = (-adiff_final LT -10))
;; chi^2 CR rejection. The CRs have to be positive, which keeps from 
;; us from masking persistence which can be negative. 
X2_adiff = adiff_cr*sqrt(ivar_adiff)
;cr_mask1 = lonarr(nx, ny)
;cr_mask2 = lonarr(nx, ny)
;FOR iorder = 0L, norders-1L DO BEGIN
;    ipix = WHERE(ordermask EQ order_vec[iorder])
;    djs_iterstat, X2_adiff[ipix], mean = mean, sigma = sigma, sigrej = 2.0
;    cr_mask1[ipix] = ((X2_adiff[ipix] GT (mean + CR_thresh*sigma)) AND $
;                      (a1[ipix] GT CR_VAL))
;    cr_mask2[ipix] = ((X2_adiff[ipix] LT (mean - CR_thresh*sigma)) AND $
;                      (a2[ipix] GT CR_VAL))
;ENDFOR
ipix = WHERE(ordermask GT 0)
djs_iterstat, X2_adiff[ipix], median = median, sigma = sigma, sigrej = 2.0
cr_mask1 = cr_sharp1 OR (X2_adiff GT  median + CR_thresh*sigma AND a1 GT CR_VAL)
cr_mask2 = cr_sharp2 OR (X2_adiff LT  median - CR_thresh*sigma AND a2 GT CR_VAL)

cr_mask1 = long(cr_mask1 EQ 0)
cr_mask2 = long(cr_mask2 EQ 0)

RETURN
END
