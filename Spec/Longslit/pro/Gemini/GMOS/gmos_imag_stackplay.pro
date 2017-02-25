

path = '/Users/joe/GMOS_redux/'
cd, path

resamp_path = path + 'resamp/'
rawpath = path + 'Raw/'
irafpath = path + 'IRAF/'

imgfiles = resamp_path + ['mgproc_imag-S20130905S0044.0001.resamp.fits' $
                          , 'mgproc_imag-S20130905S0046.0001.resamp.fits' $
                          , 'mgproc_imag-S20130905S0048.0001.resamp.fits' $
                          , 'mgproc_imag-S20130905S0051.0001.resamp.fits' $
                          , 'mgproc_imag-S20130905S0053.0001.resamp.fits' $
                          , 'mgproc_imag-S20130905S0055.0001.resamp.fits']

nimgs = n_elements(imgfiles)

if (NOT keyword_set(sigrej)) then begin
   if (nimgs LE 2) then sigrej = 1.0 $ ; Irrelevant for only 1 or 2 files
   else if (nimgs EQ 3) then sigrej = 1.1 $
   else if (nimgs EQ 4) then sigrej = 1.3 $
   else if (nimgs EQ 5) then sigrej = 1.6 $
   else if (nimgs EQ 6) then sigrej = 1.9 $
   else sigrej = 2.0
endif
if (NOT keyword_set(maxiter)) then maxiter = 3


obj = mrdfits(imgfiles[0], 0)
size = size(obj, /dim)
nx = size[0]
ny = size[1]
img_stack = fltarr(nx, ny, nimgs)
msk_stack = fltarr(nx, ny, nimgs)
var_stack = fltarr(nx, ny, nimgs)
weights   = fltarr(nx, ny, nimgs)

FOR ii = 0L, nimgs-1L DO BEGIN
   img_stack[*, *, ii] = mrdfits(imgfiles[ii], 0)
   ;; for this mask 0=good, bad=1 
   msk_stack[*, *, ii] = (img_stack[*, *, ii] LT 1e-5) 
   ;; kludge until I get a mask. Note here I reversed the masking convention 
   ;;var_stack[*, *,ii] = ???
   weights[*, *, ii] = 1.0 ;; equal weightin for testing, later this is (S/N)^2
ENDFOR

IF nimgs GT 1 THEN BEGIN
;; Average the images with rejection. Mask convention: 0 = good, 1=bad
   imgfinal_avs = djs_avsigclip(img_stack, 3, sigrej = sigrej $
                                , maxiter = maxiter, inmask = msk_stack $
                                , outmask = outmask_img)
;; Combine the variances
   nused = total(outmask_img EQ 0, 3)
   varfinal_avs = total(var_stack*(outmask_img EQ 0), 3)/(nused^2 + (nused EQ 0))
   maskfinal = (total(outmask_img, 3) NE nimgs)
; Optimally combine the images
   weights = weights*float(outmask_img EQ 0)
   wght_sum = total(weights, 3)
   imgfinal = total(weights*img_stack, 3)/(wght_sum + (wght_sum EQ 0.0))
   varfinal = total(weights^2*var_stack, 3)/(wght_sum + (wght_sum EQ 0.0))^2
ENDIF ELSE BEGIN
   maskfinal = (msk_stack EQ 0)
   imgfinal = img_stack
   varfinal = var_stack
ENDELSE


END


