PRO LONG_COADD2D, scifiles, outfile = outfile, weights = weights1 $
                  , imgfinal = imgfinal, subfinal = subfinal  $
                  , ivarfinal = ivarfinal, maskfinal = maskfinal $
                  , IREF = IREF, YMULT = YMULT, WAVEIMG = WAVEIMG


  IF n_elements(iref) EQ 0 THEN iref = 0
  nimgs = n_elements(scifiles)
  scihdr = headfits(scifiles[iref])
  objstruct = xmrdfits(scifiles[iref], 5)

IF NOT KEYWORD_SET(SIGREJ) THEN SIGREJ = 3.0
;if (NOT keyword_set(sigrej)) then begin
;   if (nimgs LE 2) then sigrej = 1.0 $ ; Irrelevant for only 1 or 2 files
;   else if (nimgs EQ 3) then sigrej = 1.1 $
;   else if (nimgs EQ 4) then sigrej = 1.3 $
;   else if (nimgs EQ 5) then sigrej = 1.6 $
;   else if (nimgs EQ 6) then sigrej = 1.9 $
;   else sigrej = 2.0
;endif
if (NOT keyword_set(maxiter)) then maxiter = 3

IF KEYWORD_SET(WEIGHTS1) THEN BEGIN
   weights = weights1
   splog, 'Using input weights'
   forprint, replicate('IMG:', nimgs), lindgen(nimgs), weights1, textout = 2
ENDIF ELSE BEGIN 
   WEIGHTS = fltarr(nimgs) + 1.0
   splog, 'No weights input, using uniform weighting'
ENDELSE
FOR ii = 0L, nimgs-1L DO BEGIN
   img  = mrdfits(scifiles[ii], 0)
   ivar = mrdfits(scifiles[ii], 1)
   sky  = mrdfits(scifiles[ii], 2)
   msk  = mrdfits(scifiles[ii], 4)
   IF ii EQ 0 THEN BEGIN
      size = size(img, /dim)
      nx = size[0]
      ny = size[1]
      img_stack = fltarr(nx, ny, nimgs)
      msk_stack = fltarr(nx, ny, nimgs)
      var_stack = fltarr(nx, ny, nimgs)
      sub_stack = fltarr(nx, ny, nimgs)
      wgt_stack = fltarr(nx, ny, nimgs)
   ENDIF
   ;; ???? In the future implement re-scaling code here. Re-scale each slit 
   ;; ???? by sky level, or perhaps the entire image. 
   ;; ???? Currently re-scaling by counts in the brightest object
   yml_now = replicate(1.0, nx) # ymult[*, ii]
   img_stack[*, *, ii] = img*yml_now
   sub_stack[*, *, ii] = (img - sky)*yml_now
   var_stack[*, *, ii] = float(ivar GT 0.0)*yml_now^2/(ivar + (ivar LE 0.0))
   ;; for this mask 0=good, bad=1 
   msk_stack[*, *, ii] = (msk EQ 0)
   ;; equal weighting for testing, later this is (S/N)^2 of the
   ;; highest S/N traces.
   wgt_stack[*, *, ii] = weights[ii]
ENDFOR

IF nimgs GT 1 THEN BEGIN
;; Average the images with rejection. Mask convention: 0 = good, 1=ba
   imgfinal_avs = djs_avsigclip(img_stack, 3, sigrej = sigrej $
                                , maxiter = maxiter, inmask = msk_stack $
                                , outmask = outmask_img)
   subfinal_avs = djs_avsigclip(sub_stack, 3, sigrej = sigrej $
                                , maxiter = maxiter, inmask = msk_stack $
                                , outmask = outmask_sub)
;; Combine the variances
   nused = total(outmask_img EQ 0, 3)
   maskfinal = (total(outmask_img, 3) NE nimgs)
;; Optimally combine the images. We use the img mask (not the sub mask)
   wgt_stack = wgt_stack*float(outmask_img EQ 0)
   wgt_sum = total(wgt_stack, 3)
   imgfinal = total(wgt_stack*img_stack, 3)/(wgt_sum + (wgt_sum EQ 0.0))
   subfinal = total(wgt_stack*sub_stack, 3)/(wgt_sum + (wgt_sum EQ 0.0))
   varfinal = total(wgt_stack^2*var_stack, 3)/(wgt_sum + (wgt_sum EQ 0.0))^2
ENDIF ELSE BEGIN
   maskfinal = (msk_stack EQ 0)
   imgfinal = img_stack
   subfinal = sub_stack
   varfinal = var_stack
ENDELSE


ivarfinal = float(maskfinal GT 0)*float(varfinal GT 0)/(varfinal + (varfinal EQ 0.0))

IF KEYWORD_SET(OUTFILE) THEN BEGIN
   sxaddpar, scihdr, 'NEXP', nimgs, ' number of exposures combined'
   sxaddpar, scihdr, 'COADD_2D', 1L, ' If true, this was a 2d coadd'
   ;;IF KEYWORD_SET(WAVEIMG) THEN naxis = 5 ELSE naxis = 4
   ;naxis = 4
   ;sxaddpar, scihdr, 'NAXIS', naxis
   mwrfits, float(imgfinal), outfile, scihdr, /create
   mwrfits, float(ivarfinal), outfile
   mwrfits, float(subfinal), outfile
   mwrfits, float(maskfinal), outfile
   ;;stop
   ;;mwrfits, objstruct, outfile
   ;;IF KEYWORD_SET(WAVEIMG) THEN mwrfits, float(waveimg), outfile
ENDIF
   
END
