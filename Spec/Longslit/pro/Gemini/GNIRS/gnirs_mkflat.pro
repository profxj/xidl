PRO gnirs_mkflat, bluefiles, redfiles, order5files, order6files, objfiles $
                  , VERBOSE = VERBOSE, tset_slits = tset_slits $
                  , flat = final_median, CONV = CONV, CHK = CHK $
                  , SHIFT8 = SHIFT8_1

  main_hdr = headfits(bluefiles[0], exten = 0, /silent)
   IF KEYWORD_SET(SHIFT8_1) THEN SHIFT8 = SHIFT8_1 $
   ELSE SHIFT8 = 0
   chk = 1
   sz = size(bluefiles)
   nblue = sz[1]
   temp = xmrdfits(bluefiles[0], 1, hdr, /silent)
   sz = size(temp)
   nx = sz[1]
   ny = sz[2]
   blue_stack = fltarr(nblue, nx, ny)

   colors=getcolor(/Load)

   print, "gnirs_mkflat: reading in BLUE flats..."
   for i=0, nblue-1 do begin
       blue_stack[i,*,*] = xmrdfits(bluefiles[i], 1, /silent)
   endfor
   
;;;;;; Read in the Red flats

   
   sz = size(redfiles)
   nred = sz[1]
   temp = xmrdfits(redfiles[0], 1, hdr, /silent)
   sz = size(temp)
   if ((nx NE sz[1]) OR (ny NE sz[2])) then begin
       print, "Image sizes do not match!!! returning..."
       stop
   endif
   red_stack = fltarr(nred, nx, ny)

   print, "gnirs_mkflat: reading in RED flats..."
   for i=0, nred-1 do begin
       red_stack[i,*,*] = xmrdfits(redfiles[i], 1, /silent)
   endfor

;;;;;; Read in the order 3 & 4 files
   
   sz = size(order5files)
   n3 = sz[1]
   temp = xmrdfits(order5files[0], 1, hdr, /silent)
   sz = size(temp)
   if ((nx NE sz[1]) OR (ny NE sz[2])) then begin
       print, "Image sizes do not match!!! returning..."
       stop
   endif
   order5_stack = fltarr(n3, nx, ny)

   print, "gnirs_mkflat: reading in Order 5 flats..."
   for i = 0, n3-1 do begin
       order5_stack[i,*,*] = xmrdfits(order5files[i], 1, /silent)
   endfor
   
   sz = size(order6files)
   n4 = sz[1]
   temp = xmrdfits(order6files[0], 1, hdr, /silent)
   sz = size(temp)
   if ((nx NE sz[1]) OR (ny NE sz[2])) then begin
       print, "Image sizes do not match!!! returning..."
       stop
   endif
   order6_stack = fltarr(n4, nx, ny)
   
   print, "gnirs_mkflat: reading in Order 6 flats..."
   for i = 0, n4-1 do begin
       order6_stack[i, *, *] = xmrdfits(order6files[i], 1, /silent)
   endfor
   
;  Read in the object file to get a wavelenght map

;   scihdr = headfits(objfile)
;   slit_hdr = strtrim(sxpar(hdr, 'SLIT'))
;   slit_str = strsplit(slit_hdr, 'arcsec', /extr)
;   slit = double(slit_str)
;   pkwdth = slit/0.15D          ; GNIRS plate scale is 0.15"
;   TOLER = pkwdth/2.0D
;   FWHM = pkwdth
;;;;;;;;;;;;;; Shift to find order edges
   print, "gnirs_mkflat: Finding orders..."
   bluemedian = djs_avsigclip(blue_stack, 1)
   IF conv EQ '2005' THEN BEGIN
       redmedian = djs_avsigclip(red_stack, 1)
       order5_median = djs_avsigclip(order5_stack, 1)
       order6_median = djs_avsigclip(order6_stack, 1)
       bluemedian    = djs_avsigclip(blue_stack, 1)
       ;; kludge this for now. Gemini took stupid calibrations
       ;; which do not cover order 8. This is so we can trace something.
       old_tsetfile = getenv('LONGSLIT_DIR') $
         + '/calib/flats/GNIRS/gnirs_archive_orders_2005.fits'
       tset_slits_old = mrdfits(old_tsetfile, 1)
       ordermask1 = long_slits2mask(tset_slits_old, xshift = -30)
       ordermask2 = long_slits2mask(tset_slits_old, xshift = 25)
       ordermask3 = long_slits2mask(tset_slits_old, xshift = -15)
       ordermask1[WHERE(ordermask1 GT 0)] = $
         ordermask1[WHERE(ordermask1 GT 0)] + 2L
       ordermask2[WHERE(ordermask2 GT 0)] = $
         ordermask2[WHERE(ordermask2 GT 0)] + 2L
       ordermask3[WHERE(ordermask3 GT 0)] = $
         ordermask3[WHERE(ordermask3 GT 0)] + 2L
       bluemedian = bluemedian*(ordermask1 EQ 8 OR $
                                ordermask2 EQ 8 OR $
                                ordermask3 EQ 8)
       bluemedian = redmedian >  order5_median >  order6_median > bluemedian
   ENDIF
;;;;;;;;;;  Now fit the vertical sense of the sawtooth image
   
;;;;;;;;;;;;;;;; Fit the left hand edge   
   
   buff = 15.0
   tset_slits_init = gnirs_traceorders(bluemedian, main_hdr)
   
   ;; Kludge to adjust last order position in 2004 data
   IF conv EQ '2005' THEN $
      tset_slits_init.coeff[0, 5] = tset_slits_init.coeff[0, 5] + shift8 ;; - 8.0

   dimt = size(tset_slits_init[0].COEFF, /dim)
   norders = dimt[1]

   ;; Expand the slit boundary for making the flat and then we will
   ;; tweak them later
   tset_slits_init[0].coeff[0, *] = tset_slits_init[0].coeff[0, *] -2.0  
   tset_slits_init[1].coeff[0, *] = tset_slits_init[1].coeff[0, *] +2.0  
   ximg = long_slits2x(tset_slits_init, edgmask = edgmask, TOL_EDG = 6.0D)

   ordermask = long_slits2mask(tset_slits_init)
   ordermask[WHERE(ordermask GT 0)] = ordermask[WHERE(ordermask GT 0)] + 2L

   
   
;  Generate a mask giving position on the slit
;  edge mask is agressive to avoid artifacts in flat
   order_vec = [3, 4, 5, 6, 7, 8]
;;  Trace the object file to get a wavelength map.
   gnirs_proc, objfiles[0], a1, ivar_a1, hdr = scihdr, /pattern
   gnirs_proc, objfiles[1], b1, ivar_b1, /pattern
   gnirs_proc, objfiles[2], b2, ivar_b2, /pattern
   gnirs_proc, objfiles[3], a2, ivar_a2, /pattern
   img_stack = fltarr(nx*ny, 4)
   maskstack = lonarr(nx*ny, 4)
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
   delvarx, img_stack
   delvarx, maskstack
   waveimg = gnirs_waveimg(wpix_img, scihdr, tset_slits_init, piximg = piximg $
                           , QAFILE = 'flat.ps')
   spawn, 'rm -f flat*'
   ;;xshift = long_xcorr_slits(wpix_img, tset_slits_init, /shift)
;; Generate the median flatfield images, both blue and red
   final_flat = fltarr(nx, ny)
   
;;;;;;;;;;  normalize out the shape of the lamp spectrum for each image
   
   print, "gnirs_mkflat: normalizing by lamp spectrum..."
;   yjunk = replicate(1,nx) # findgen(ny) 

   CASE CONV OF
       '2004': BEGIN
           redorders = [3, 4]
           order5orders = [5]
           order6orders = [6]
           blueorders = [7, 8]
       END
       '2005': BEGIN
           redorders = [3, 4, 5]
           order5orders = [6]
           order6orders = [7]
           blueorders = [8]
       END
       '2006_old': BEGIN
           redorders    = [ 3, 4 ]
           order5orders = [ 5 ]
           order6orders = [ 6 ]
           blueorders   = [ 7, 8]
       END
       '2007': BEGIN
           redorders    = [ 3 ]
           order5orders = [ 4, 5 ]
           order6orders = [ 6 ]
           blueorders   = [ 7, 8] 
       END
   ENDCASE
   
;      fit each member of the blue stack
   FOR k = 0L, nblue -1L DO BEGIN
       thisimg = reform(blue_stack[k, *, *], nx*ny)
       FOR iorder = 0L, n_elements(blueorders)-1L DO BEGIN
           order = blueorders[iorder]
           this_order = WHERE(ordermask EQ order)
           gd = where(ordermask EQ order AND $
                      ximg GE 0.1 AND ximg LE 0.9, npts)
           pix = piximg[gd]
           fxarr = thisimg[gd]
           fxarr = fxarr[sort(pix)]
           pix = pix[sort(pix)]
           
           med_fx = djs_median(fxarr, width = 20)
           fitpix = WHERE(fxarr GE med_fx-3*sqrt(med_fx) AND $
                          fxarr LE med_fx+3*sqrt(med_fx), nfit)
           
           sset = bspline_iterfit(pix[fitpix], fxarr[fitpix] $
                                  , nbkpts = ny/10 $
                                  , upper = 2, lower = 2, nord = 3 $
                                  , yfit = bfit)
           ymodel = bspline_valu(piximg[this_order], sset)
           
           IF KEYWORD_SET(CHK) THEN BEGIN
               plot, pix[fitpix], fxarr[fitpix], psym = 3 $ 
                     , xrange = [min(pix[fitpix]), max(pix[fitpix])] $
                     , xstyle = 1
               oplot, pix[fitpix], bfit, color = djs_icolor('red') ;;, 245)
               wait,1.0
           ENDIF
           thisimg[this_order] = thisimg[this_order]/ymodel
       ENDFOR
       blue_stack[k, *, *] = reform(thisimg, nx, ny)
   ENDFOR
   
;      Fit each member of the red stack       
   FOR k = 0L, nred -1L DO BEGIN
       thisimg = reform(red_stack[k, *, *], nx*ny)
       FOR iorder = 0L, n_elements(redorders)-1L DO BEGIN
           order = redorders[iorder]
           this_order = WHERE(ordermask EQ order)
           gd = where(ordermask EQ order  $
                      AND ximg GE 0.1 AND ximg LE 0.9, npts)
           pix = piximg[gd]
           fxarr = thisimg[gd]
           fxarr = fxarr[sort(pix)]
           pix = pix[sort(pix)]
           
           med_fx = djs_median(fxarr, width = 20)
           fitpix = WHERE(fxarr GE med_fx-3*sqrt(med_fx) AND $
                          fxarr LE med_fx+3*sqrt(med_fx), nfit)
           
           sset = bspline_iterfit(pix[fitpix], fxarr[fitpix] $
                                  , nbkpts = ny/10 $
                                  , upper = 2, lower = 2, nord = 3 $
                                  , yfit = bfit)
           ymodel = bspline_valu(piximg[this_order], sset)
           
           IF KEYWORD_SET(CHK) THEN BEGIN
               plot, pix[fitpix], fxarr[fitpix], psym = 3 $ 
                     , xrange = [min(pix[fitpix]), max(pix[fitpix])] $
                     , xstyle = 1
               oplot, pix[fitpix], bfit, color = djs_icolor('red');;, 245)
               wait,1.0
           ENDIF
           thisimg[this_order] = thisimg[this_order]/ymodel
       ENDFOR
       red_stack[k, *, *] = reform(thisimg, nx, ny)
   ENDFOR
;      Fit each member of the order5 stack
   FOR k = 0L, n3-1L DO BEGIN
       thisimg = reform(order5_stack[k, *, *], nx*ny)
       FOR iorder = 0L, n_elements(order5orders)-1L DO BEGIN
           order = order5orders[iorder]
           this_order = WHERE(ordermask EQ order)
           gd = where(ordermask EQ order  $
                      AND ximg GE 0.1 AND ximg LE 0.9, npts)
           pix = piximg[gd]
           fxarr = thisimg[gd]
           fxarr = fxarr[sort(pix)]
           pix = pix[sort(pix)]
           
           med_fx = djs_median(fxarr, width = 20)
           fitpix = WHERE(fxarr GE med_fx-3*sqrt(med_fx) AND $
                          fxarr LE med_fx+3*sqrt(med_fx), nfit)
           
           sset = bspline_iterfit(pix[fitpix], fxarr[fitpix] $
                                  , nbkpts = ny/10 $
                                  , upper = 2, lower = 2, nord = 3 $
                                  , yfit = bfit)
           ymodel = bspline_valu(piximg[this_order], sset)
           
           IF KEYWORD_SET(CHK) THEN BEGIN
               plot, pix[fitpix], fxarr[fitpix], psym = 3 $ 
                     , xrange = [min(pix[fitpix]), max(pix[fitpix])] $
                     , xstyle = 1
               oplot, pix[fitpix], bfit, color = djs_icolor('red');, 245)
               wait,1.0
           ENDIF
           thisimg[this_order] = thisimg[this_order]/ymodel
       ENDFOR
       order5_stack[k, *, *] = reform(thisimg, nx, ny)
   ENDFOR
   FOR k = 0L, n4-1L  DO BEGIN
       thisimg = reform(order6_stack[k, *, *], nx*ny)
       FOR iorder = 0L, n_elements(order6orders)-1L DO BEGIN
           order = order6orders[iorder]
           this_order = WHERE(ordermask EQ order)
           gd = where(ordermask EQ order  $
                      AND ximg GE 0.1 AND ximg LE 0.9, npts)
           pix = piximg[gd]
           fxarr = thisimg[gd]
           fxarr = fxarr[sort(pix)]
           pix = pix[sort(pix)]
           
           med_fx = djs_median(fxarr, width = 20)
           fitpix = WHERE(fxarr GE med_fx-3*sqrt(med_fx) AND $
                          fxarr LE med_fx+3*sqrt(med_fx), nfit)
           
           sset = bspline_iterfit(pix[fitpix], fxarr[fitpix] $
                                  , nbkpts = ny/10 $
                                  , upper = 2, lower = 2, nord = 3 $
                                  , yfit = bfit)
           ymodel = bspline_valu(piximg[this_order], sset)
           IF KEYWORD_SET(CHK) THEN BEGIN
               plot, pix[fitpix], fxarr[fitpix], psym = 3 $ 
                     , xrange = [min(pix[fitpix]), max(pix[fitpix])] $
                     , xstyle = 1
               oplot, pix[fitpix], bfit, color = djs_icolor('red');;, 245)
               wait, 1.0
           ENDIF
           thisimg[this_order] = thisimg[this_order]/ymodel
       ENDFOR
       order6_stack[k, *, *] = reform(thisimg, nx, ny)
   ENDFOR
   
   if (nblue LE 2) then sigrej_blue = 1.0 $ 
   else if (nblue EQ 3) then sigrej_blue = 1.1 $
   else if (nblue EQ 4) then sigrej_blue = 1.3 $
   else if (nblue EQ 5) then sigrej_blue = 1.6 $
   else if (nblue EQ 6) then sigrej_blue = 1.9 $
   else sigrej_blue = 2.0
   
   if (nred LE 2) then sigrej_red = 1.0 $ 
   else if (nred EQ 3) then sigrej_red = 1.1 $
   else if (nred EQ 4) then sigrej_red = 1.3 $
   else if (nred EQ 5) then sigrej_red = 1.6 $
   else if (nred EQ 6) then sigrej_red = 1.9 $
   else sigrej_red = 2.0
   
   if (n3 LE 2) then sigrej_3 = 1.0 $ 
   else if (n3 EQ 3) then sigrej_3 = 1.1 $
   else if (n3 EQ 4) then sigrej_3 = 1.3 $
   else if (n3 EQ 5) then sigrej_3 = 1.6 $
   else if (n3 EQ 6) then sigrej_3 = 1.9 $
   else sigrej_3 = 2.0
   
   if (n4 LE 2) then sigrej_4 = 1.0 $ 
   else if (n4 EQ 4) then sigrej_4 = 1.1 $
   else if (n4 EQ 4) then sigrej_4 = 1.4 $
   else if (n4 EQ 5) then sigrej_4 = 1.6 $
   else if (n4 EQ 6) then sigrej_4 = 1.9 $
   else sigrej_4 = 2.0
   
   
   blue_median = djs_avsigclip(blue_stack, 1, sigrej = sigrej_blue $
                               , inmask = (ordermask EQ 0) $
                               , maxiter = maxiter)
   red_median = djs_avsigclip(red_stack, 1, sigrej = sigrej_red $
                              , inmask = (ordermask EQ 0) $
                              , maxiter = maxiter)
   order5_median = djs_avsigclip(order5_stack, 1, sigrej = sigrej_3 $
                                 , inmask = (ordermask EQ 0) $
                                 , maxiter = maxiter)
   order6_median = djs_avsigclip(order6_stack, 1, sigrej = sigrej_4 $
                                 , inmask = (ordermask EQ 0) $
                                 , maxiter = maxiter)
   
   
;     combine into one flat
   FOR iorder = 0L, norders-1L do begin
       order = order_vec[iorder]
       this_order = WHERE(ordermask EQ order)
       if (max(order EQ redorders) EQ 1) then begin
          final_flat[this_order] = red_median[this_order]
       endif 
       
       if (max(order EQ blueorders) EQ 1) then begin
          final_flat[this_order] = blue_median[this_order]
       endif
       
       if (max(order EQ order5orders) EQ 1) then begin
          final_flat[this_order] = order5_median[this_order]
       endif
       
       if (max(order EQ order6orders) EQ 1) then begin
          final_flat[this_order] = order6_median[this_order]
       endif
    ENDFOR

   ;; Now tweak slit boundaries
   chk = 1
   tset_slits = long_slittweak(final_flat, tset_slits_init, CHK = chk)
   ordermask = long_slits2mask(tset_slits)
   ordermask[WHERE(ordermask GT 0)] = ordermask[WHERE(ordermask GT 0)] + 2L
   final_median = fltarr(nx, ny) + 1.0

   ;     combine into one flat
   FOR iorder = 0L, norders-1L do begin
       order = order_vec[iorder]
       this_order = WHERE(ordermask EQ order)
       if (max(order EQ redorders) EQ 1) then begin
          final_median[this_order] = red_median[this_order]
       endif 
       
       if (max(order EQ blueorders) EQ 1) then begin
          final_median[this_order] = blue_median[this_order]
       endif
       
       if (max(order EQ order5orders) EQ 1) then begin
          final_median[this_order] = order5_median[this_order]
       endif
       
       if (max(order EQ order6orders) EQ 1) then begin
          final_median[this_order] = order6_median[this_order]
       endif
    ENDFOR
   
   yimg = replicate(1.0, nx) # findgen(ny)
   
   IF conv EQ '2005' THEN BEGIN
       ;; Since Gemini took stupid calibraitons we only flat field
       ;; orders 3-5
       unitlogic = (ordermask EQ 3 AND waveimg LT 19300.0D) OR $
         (ordermask EQ 4 AND waveimg GT 18050.0D) OR $
         (ordermask EQ 4 AND yimg GT 1000) OR $
         (ordermask EQ 5 AND waveimg LT 14160.0D) OR $
         (ordermask EQ 6) OR (ordermask EQ 7) OR (ordermask EQ 8) OR $
         edgmask EQ 1
   ENDIF ELSE BEGIN
       unitlogic = (ordermask EQ 3 AND waveimg LT 19300.0D) OR $
         (ordermask EQ 4 AND waveimg GT 18050.0D) OR $
         (ordermask EQ 4 AND yimg GT 1000) OR $
         (ordermask EQ 5 AND waveimg GT 13290.0D AND waveimg LT 14160.0D) OR $
         (ordermask EQ 5 AND yimg GT 875) OR $
         (ordermask EQ 6 AND yimg GT 800) OR $
         (ordermask EQ 7 AND yimg GT 740) OR $
         (ordermask EQ 8) OR $
         edgmask EQ 1
   ENDELSE
   unitinds = WHERE(unitlogic, nunit)
   ;;     (ordermask EQ 8 AND yimg GT 590) OR $
   IF nunit GT 0 THEN final_median[unitinds] = 1.0
   
;   splog, 'Writing output file'
;   mwrfits, final_median, outfile, hdr, /create
;   mwrfits, tset_slits, outfile
       
;    tmp = { flatstruct }
;    tmp.flat = final_median
;    tmp.order_img = order_image
;    tmp.norders = norders
;    tmp.nx = nx
;    tmp.ny = ny
;    tmp.leftedge = leftedge
;    tmp.rightedge = rightedge
   
;   return, tmp
   
   print, ""
   print, "gnirs_mkflat: all done!"
   print, ""

   RETURN   
END
