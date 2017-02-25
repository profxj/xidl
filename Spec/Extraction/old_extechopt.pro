;+ 
; NAME:
; x_extechopt
;     Version 1.1
;
; PURPOSE:
;   Extract 1D spectra from the 2D images.  For each order, a boxcar
;   and an optimal extraction is performed.  For the latter, an object
;   profile is derived and both the object flux and sky are fit throughout
;   the order (i.e. not row by row).  The main driver is
;   x_extechopt.  For the optimal extraction, the data is extracted
;   to a specific set of vacuum wavelengths, chosen to be the same for
;   every spectrum to facilitate coadding without rebinning.  Here are
;   the steps in detail:
;
;    1.  Perform a boxcar extraction using extract_box
;    2.  Estimate the SNR per order from the boxcar extraction
;      -- LOOP ON ORDERS IN DECREASING SNR --
;    3.  Fit the boxcar extraction with a bspline
;    4.  Calculate the object profile
;          a. bspline_iterfit the flux vs position on slit
;          b. Force the profile to be positive everywhere and have a
;          sensible FWHM
;    5.  Fit the order using the profile and sky (bspline_extract)
;
; CALLING SEQUENCE:
;   
; x_extechopt, img, skysub, ivar, ordr_str, obj_str, velpix, radius=
;
; INPUTS:
;  img      -- Data image
;  skysub   -- Sky subtraction image from x_echskysub
;  ivar     -- Inverse variance image
;  ordr_str -- Order structure describing the echelle footprint
;  velpix   -- Size of pixel in velocity units (km/s)
;
; RETURNS:
;
; OUTPUTS:
;  obj_str  -- Structure containing the trace and extracted data
;
; OPTIONAL KEYWORDS:
;   /BOXONLY  -- Only do boxcar extraction
;   SKYFIL=  -- Filename containing the sky spectrum
;   IMG_ARC=  -- Wavelength image
;   BASE_APER= -- Aperture for boxcar [default: 0.75, 0.75]
;   /OCHK   - Plot the extracted flux (optimal) for each order
;   /RESCHK - Check the residuals of the 2D, fully extracted image
;   /DEBUG  - Stop within extraction routine to check stuff
;   HIGHSNR - Value of SNR^2 of the data for a given order which when
;             exceeded mike_box uses an additional parameter for the
;             profile shape.  (Default:  500 corresponding to SNR=22)
;             Lowering this parameter may improve extraction.
;   ORDRS   - Input array of physical order numbers to extract
;   /EXTENBOX   - Allows the boxcar aperture to be larger than the 
;                 slit length defined by the trace flats (only
;                 recommended for very bright stars)
;   MIN_CUT - Minimum length for cutoff aper in fractional slit length
;             [default: 0.].  For
;             longer slits and good seeing, it may be useful to set
;             this to a number like 0.5 or larger
;   FIN_MASK -- Image of masked pixels (mainly cosmic rays)
;   MSKTRIM= -- Value to trim edge of slit by [default: 0.5]
;   /SKIPSKYSUB  -- Skip sky subtraction
;   /NOCRMASK - ignore the CR flagging 
;
; Optional OUTPUTS:
;  FIN_TRC=  -- Final trace of spectrum
;  MODEL_OBJ= -- Model of the object flux
;  MODEL_SKY= -- Model of the sky
;  MODEL_PROF= -- Model of the object profile
;
; COMMENTS:
;  The program extracts the orders in order of decreasing SNR.  If the
;  SNR is lower than lowsnr (default: 2.49) then the optimal
;  extraction is performed using the profile parameters from the
;  previous order(s).
;
; EXAMPLES:
;   x_extechopt
;
; PROCEDURES/FUNCTIONS CALLED:
;  extract_boxcar
;  x_smoothmask
;  bspline_extract
;
; REVISION HISTORY:
;   26-Aug-2003 Written by SMB
;   Feb-2005 Ported to XIDL by JXP
;   2009-Jan-03 Moustakas - added NOCRMASK
;-
;------------------------------------------------------------------------------

function x_smoothmask, outmask, indx, ncol
   
     if N_PARAMS() LT 3 then message, 'x_smoothmask(outmask, indx, ncol)' 
     x = indx mod ncol
     y = indx /  ncol

     nx = max(x) - min(x) + 5L
     ny = max(y) - min(y) + 5L
     new_indx = (y-min(y)+2) *nx + x - min(x) + 2 
     
     temp_arr  = fltarr(nx,ny)
     temp_arr[new_indx] = 1.0 * (outmask EQ 0)
    
     kernel = [[0,0,1,0,0],[0,1,1,1,0],[1,1,1,1,1],[0,1,1,1,0],[0,0,1,0,0]] 
     smooth_arr = convol(temp_arr,kernel)
     
     return_mask = smooth_arr[new_indx] EQ 0
   
return, return_mask
end 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_extechopt, img, skysub, ivar, ordr_str, obj_str, velpix, $
                 base_aper=base_aper, chk=chk, img_arc=img_arc, $
                 OCHK=ochk, RESCHK=reschk, HIGHSNR=highsnr, $
                 ORDRS=ordrs, skyfil=skyfil, BOXONLY=boxonly, $
                 helio=helio, DEBUG=debug, OBJ_NAME=obj_name, $
                 FIN_TRC=fin_trc, ORDERMASK = ordermask, $
                 MODEL_OBJ=model_obj, MODEL_SKY=model_sky, SLIT_LEN=slit_len, $
                 MODEL_PROF=model_prof, MSKTRIM=msktrim, EXTENBOX=extenbox,$
                 SKIPSKYSUB=skipskysub, MIN_CUT=min_cut, $
                 EXTRACT_mask=extract_mask, TST_IVAR=tst_ivar, nocrmask=nocrmask ; jm09jan03nyu

  if  N_params() LT 5  then begin 
      print,'Syntax - ' + $
        'x_extechopt, img, skysub, ivar, ordr_str, obj_str, velpix', $
        '/DEBUG, /CHK, /OCHK, /RESCHK, /READSKY, ORDRS=,'
      print, '          HIGHSNR=, MIN_CUT=, /SKIPSKYSUB [v1.1]'
      return
  endif 
  

  if not keyword_set(SLIT_LEN) then slit_len = 5.0 ; slit length in arcseconds
  if NOT keyword_set(MSKTRIM) then msktrim = 0.5
  if NOT keyword_set(MIN_CUT) then min_cut = 0.
  if NOT keyword_set(HIGHSNR) then highsnr = 500.
  if NOT keyword_set(skyfil) then skyfil=''
  if N_ELEMENTS(FIN_TRC) EQ 0 then fin_trc = 2

  inverse_bin = (size(skysub))[1]/2048.0
  if NOT keyword_set(base_aper) then base_aper=[0.75,0.75]
  if NOT keyword_set(radius) then begin
      if keyword_set(mean(obj_str.aper)) then radius=mean(obj_str.aper) $
      else radius = 10.*inverse_bin
  endif
  
  if NOT keyword_set(obj_name) then obj_name = ' '
  if not keyword_set(BSPLINE_ORD) then bspline_ord = 4L 
;  if not keyword_set(BSPLINE_ORD) then bspline_ord = 4L 
  
  ;; first is ultra-easy boxcar extraction..
  sz = size(skysub, /dimen)
  
  ncol = sz[0]
  nrow = sz[1]
  
  badpix = where((finite(skysub) EQ 0) OR (finite(ivar) EQ 0) $
                 OR ivar LE 0., nbadpix)
  if nbadpix GT 0 then begin
;       print, "NaN's in this image!", nbadpix
      skysub[badpix] = 0
      ivar[badpix] = 0
  endif 
  var = 1/(ivar + (ivar EQ 0)) * (ivar GT 0)
  
  ordr_shift  = ordr_str
  slit_length = ordr_shift.rhedg - ordr_shift.lhedg
  half_length = slit_length/2.
  
  objcen      = obj_str.trace[0:sz[1]-1]  
  
  left_edge   = (objcen - base_aper[0] * half_length) > ordr_shift.lhedg
  right_edge  = (objcen + base_aper[1] * half_length) < ordr_shift.rhedg
  
  fx   =     extract_asymbox2(skysub, left_edge, right_edge)
  fvar  = extract_asymbox2(var, left_edge, right_edge)
  fmask = extract_asymbox2(1.0*(ivar LE 0.), left_edge, right_edge)

  fvar  = fvar*(fmask EQ 0) - (fmask GT 0)
  fi    = 1.0/(fvar + (fvar EQ 0)) * (fvar GT 0)
  median_sn1 = djs_median(fx * sqrt(fi),1)


;;;;;;;;
;          Sort from highest S/N order to lowest
;
  ;; only keep ones with (flg_anly AND 1)
  do_orders = where(obj_str.flg_anly AND 1)
  if do_orders[0] NE -1 then begin
      order_table = reverse(do_orders[sort(median_sn1[do_orders])]) 
  endif else begin
      print, 'No orders found with flg_anly AND 1'
      order_table = reverse(sort(median_sn1)) 
  endelse
  
  if NOT keyword_set(LOWSNR) then lowsnr = 2.49
  nuseable = long(total(median_sn1 GT LOWSNR))
  print, 'Found ', nuseable, ' orders with SNR greater than ', LOWSNR, $
    format='(a,i3,a,f6.2)'
  
  badpix = where(finite(fx) EQ 0)
  if badpix[0] NE -1 then begin
      fx[badpix] = 0
      fi[badpix] = 0
  endif
  
  
  ;; Shift the Edge orders
  if not keyword_set(ORDRS) then nord = n_elements(order_table) $
  else nord = n_elements(ordrs)
  
  ;; Create the order mask (img with val = order # or -order# for gaps)
  ordermask = x_ordermask(sz[0], sz[1], ordr_shift, trim=msktrim)
  
  ordr_fake = ordr_shift
  ordr_fake.lhedg = ordr_fake.lhedg+200
  ordr_fake.rhedg = ordr_fake.rhedg+200
  fakemask = x_ordermask(sz[0]+400, sz[1], ordr_fake, trim=0)
  
  ;; Create final images
  extract_mask  = ordermask*0 + 1
  model_sky     = ordermask*0.0
  model_obj     = ordermask*0.0
  model_prof    = ordermask*0.0
  chi_image     = ordermask*0.0
  medrow = sz[1]/2
  ycol = dindgen(sz[1])
  
  fwhm_here = fltarr(n_elements(obj_str)+1)
  
  ;; Loop on Orders (in decreasing SNR from boxcar)
  
  if keyword_set(chk) then begin 
      clr = getcolor()
      clr = getcolor(/load)
  endif

  FOR i=0,nord-1 do begin 
;  FOR i=5,nord-1 do begin 
      
      if not keyword_set( ORDRS ) then q = order_table[i] $ ; Standard
      else q = (where(ordr_shift.order EQ ordrs[i]))[0] ; Input order number
      ;; 
      print, 'x_extechopt: Extracting Order ', ordr_shift[q].order
      objcen = obj_str[q].fin_trc[0:sz[1]-1] 
      ;;
      ;;
      if keyword_set(img_arc) then begin
          inorder = where(ordermask EQ ordr_shift[q].order AND img_arc GT 3.) 
      endif else inorder = where(ordermask EQ ordr_shift[q].order)
      
      ;; x, y positions of the pixels in order q
      ystart  = 1.0d*(inorder / ncol)
      xstart  = 1.0d*(inorder mod ncol)  
;     aper    = obj_str[q].aper > base_aper
      aper    = base_aper
      
      y = 2*ystart/sz[1] - 1.
      slit_cen  = 0.5*(ordr_shift[q].lhedg + ordr_shift[q].rhedg)[ystart]
      slit_cen_fix  = 0.5*(ordr_shift[q].lhedg + ordr_shift[q].rhedg)
      slit_pos  = xstart - slit_cen
      
      ywave = x_qckwav(slit_pos, ystart, ordr_shift[q].arc_m, $
                       arc_slope=arc_slope)
      ;;
      ;;  Slit_frac goes from -1 to 1 
      ;;

      slit_frac = frac_order(ordr_str[q], xstart, ywave)
      slit_proj    = slitproj(ystart, ywave, arc_slope, ordr_shift[q])
      
      ;; Grab the slit profile for each pixel
      slit_profile = x_slitprofile_return(slit_frac, ystart, ordr_str[q])
      
      ;;  Let's fit vs. alog10(wavelengths) if possible
      objcen = obj_str[q].trace[0:sz[1]-1] 
      if keyword_set(img_arc) then begin
          wave_order = img_arc[inorder] 
          central_wave = extract_boxcar(img_arc, objcen, radius=0.5)
          wave_diff = wave_order - central_wave[ystart]
          zero_wave = where(wave_order LT 3.0 OR (abs(wave_diff) GT 0.0003), $
                           NCOMPLEMENT=n_nonzero)
          if zero_wave[0] NE -1 then begin
              print, 'Rejecting ', n_elements(zero_wave), $
                ' pixels which have wavelengths from neighboring orders', $
                format ='(a,i6,a)'
              ivar[inorder[zero_wave]] = 0
              wave_order[zero_wave] = 0
          endif
          flg_zero = 0
          if n_nonzero GT 0 then begin 
              raw_wave_set = bspline_iterfit(ywave, 1.0d*wave_order, $
                                             invvar=1.0d*(wave_order GT 3), $
                                             nbkpt=10, /silent, /double)
              
              
              dbl_wave_set = { fullbkpt : 1.0d*raw_wave_set.fullbkpt, $
                               bkmask   : raw_wave_set.bkmask, $
                               nord     : raw_wave_set.nord, $
                               coeff    : 1.0d*raw_wave_set.coeff, $
                               icoeff   : 1.0d*raw_wave_set.icoeff}
          endif else begin
              print, 'x_extechopt: No non-zero wavelengths.  Skipping...'
              flg_zero = 1 ;; Toss out bad order
          endelse
      endif else wave_order = ywave

      if flg_zero then continue
      
      ys    = sort(wave_order)
      
      sky_guess=0
;
;    This is quick skysubtraction, in case skysub was skipped
;
      if keyword_set(skipskysub) then begin
          maxdist = max(abs(slit_frac))
          edgepix = abs(slit_frac) GT maxdist-0.15
          usethese = where(edgepix, nuse)
          if nuse GT 50 then begin 
              sky_repl = bspline_iterfit(ywave[usethese], $
                                         skysub[inorder[usethese]], $
                                         invvar=ivar[inorder[usethese]], $
                                         everyn=4, yfit=sky_guess, /sil)
              
              sky_guess = bspline_valu(ywave, sky_repl)
              apodization= 10.
          endif
      endif
      
      For I_TRC=0, FIN_TRC do begin
          
          if (I_TRC NE 0) AND (total(obj_str[q].fin_trc) NE 0) then begin
              print, 'Iteration ', I_TRC+1, $
                ', using trace centers stored in obj_str.fin_trc', $
                format='(a,i2,a)'
              objcen = obj_str[q].fin_trc[0:sz[1]-1] 
          endif else begin
              print, 'Iteration ', I_TRC+1, $
                ', using trace centers stored in obj_str.trace', $
                format='(a,i2,a)'
              objcen = obj_str[q].trace[0:sz[1]-1]  
          endelse
          
          obj_trace =  objcen[ystart]
          obj_pos = xstart - obj_trace  ;; Offset from obj center

          ;;   
          ;;    Slit distance is position in slit relative to object trace
          ;;      in units of arcseconds (that's the 5. for MIKE)
          ;;
          slit_dist = frac_order(ordr_str[q], xstart, ywave, ocen=obj_trace)
          
          if not keyword_set(EXTENBOX) then begin
              left_edge   = (objcen - base_aper[0] * half_length[*,q]) $
                            > ordr_shift[q].lhedg[*]
              right_edge  = (objcen + base_aper[1] * half_length[*,q]) $
                            < ordr_shift[q].rhedg[*]
          endif else begin
              left_edge   = (objcen - base_aper[0] * half_length[*,q]) > 0.
              right_edge  = (objcen + base_aper[1] * half_length[*,q]) < (sz[0]-1)
          endelse
;          left_edge   = (objcen - base_aper[0] * half_length) $
;            > ordr_shift[q].lhedg
;          right_edge  = (objcen + base_aper[1] * half_length) $
;            < ordr_shift[q].rhedg
          
          fx = extract_asymbox2(skysub, left_edge, right_edge)
          fc = extract_asymbox2(skysub, left_edge+base_aper[0]/2*half_length, $
                                right_edge-base_aper[0]/2*half_length)
          fvar  = extract_asymbox2(var, left_edge, right_edge)
          fmask = extract_asymbox2(1.0*(ivar LE 0.), left_edge, right_edge)
          fvar  = fvar*(fmask EQ 0) - (fmask GT 0)
          fi    = 1.0/(fvar + (fvar EQ 0)) * (fvar GT 0)
          
          ;; Smooth the box car extraction and identify CR's
          if keyword_set(nocrmask) then begin ; jm09jan03nyu
             med_mask = byte(fx*0.0)+1B
          endif else begin
             med_filter = long(23 * sz[1] /4096.)
             med_diff = fx - median(fx, med_filter)
             djs_iterstat, med_diff, sigma=sigma_diff, mask=med_mask, sigrej=5
          endelse
          
          ;; Fit the flux avoiding the CR's
          fx_set = bspline_iterfit(ycol, fx, invvar=fi*med_mask, $
                                   everyn=2, outmask=outmask, $
                                   yfit=fx_fit, /silent)
          fx_set = bspline_iterfit(ycol, fx, invvar=fi*med_mask*outmask,$
                                   everyn=2, outmask=outmask, $
                                   yfit=fx_fit, /silent)
          
          fc_set = bspline_iterfit(ycol, fc, invvar=fi*med_mask, $
                                   everyn=2, outmask=outmask_c, $
                                   yfit=fc_fit, /silent)
          fc_set = bspline_iterfit(ycol, fc, invvar=fi*med_mask*outmask_c,$
                                   everyn=2, outmask=outmask_c, $
                                   yfit=fc_fit, /silent)
          
          if (max(abs(fx_fit)) GT 1.0e6) then begin
              print, 'WARNING: Huge counts in smooth fit to boxcar', $
                max(abs(fx_fit))
          endif
          
          obj_str[q].box_fx[0:nrow-1] = fx
          obj_str[q].box_var[0:nrow-1] = fvar*(med_mask NE 0)-(med_mask EQ 0)
          
          
          ybox = x_qckwav( objcen - slit_cen_fix, ycol, ordr_shift[q].arc_m)
          if keyword_set(img_arc) then begin
              ;; Wavelength
              obj_str[q].nrow = nrow
              obj_str[q].box_wv[0:nrow-1] = 10^bspline_valu(ybox, dbl_wave_set)
          endif else begin
              obj_str[q].nrow = nrow
              obj_str[q].box_wv[0:nrow-1] = 1.0d*ybox
          endelse
        
          ;;;;;;;;;;;;;;
          ;; Here ENDS boxcar
          if keyword_set( BOXONLY ) then continue

          ;; Estimate the S/N from the quick box car extraction
          sn2  = fx^2 * (fi*med_mask)
          sort_sn2 = sort(sn2)
          median_sn2 = (mean(sn2[sort_sn2[0.1*nrow:0.9*nrow]]))[0]
          
          final_obj_profile=0
          

          ;; Evaluate the flux at all pixels assuming the constant
          ;; wavelength contours determined from the Arc
          fx_model = bspline_valu(ywave, fx_set) 
          fc_model = bspline_valu(ywave, fc_set) 
          t = fx_model * outmask[ystart] * 2.0/slit_proj
          t_factor = fc_model * outmask_c[ystart] * 2.0/slit_proj
          
          ;; Identify the pixels with significant flux (>10) to build a profile
          profile_temp = (skysub[inorder]-sky_guess)/(t + (t EQ 0)) * (t GT 1) 
          
          profile_ivar = ivar[inorder]* t_factor^2 * (t GT 1) * $
            (profile_temp GT -4) * (profile_temp LT 20.0) 
          
          prof_mask = (profile_ivar GT 0)
;          stop
          
;
;         Some extra code to measure sigma at different positions in frame
;
;
;       ; let's find derivative of trace 
;       yspot = round((0.5+(findgen(5)-2)/6.) * nrow)
;       hw = nrow/20
;       xw = [-1.0*(reverse(findgen(aper[0])+1)),0, findgen(aper[1])+1]
;       nw = n_elements(xw)
;       model_prof[inorder] = (skysub[inorder]-sky_guess)
;       
;       for iy=0,n_elements(yspot)-1 do begin
;          yy = yspot[iy]+lindgen(2*hw+1)-hw
;          trace_slope = (objcen[yy+1,q]-objcen[yy-1,q]) * xbin/(2*ybin)
;          xt = xw # replicate(1,2*hw+1) + objcen[yy,q] ## (xw*0. + 1)
;          yt = (xw # trace_slope)*(xbin/ybin) + yy ## (xw*0. + 1)
;
;          if iy EQ 0 then begin
;            profile_rect = interpolate(model_prof, xt, yt, cubic=0.5)
;            rect_ivar = interpolate(ivar, xt, yt, cubic=0.5)
;          endif else begin
;            profile_rect = [profile_rect, $
;                      interpolate(model_prof, xt, yt, cubic=0.5)]
;            rect_ivar = [rect_ivar , interpolate(ivar, xt, yt, cubic=0.5)]
;          endelse
;       endfor
;       xspot = (findgen(5)*nw + long(aper[0])) ## replicate(1,2*hw+1)
;
;       sigma = (fltarr(5)+1.) ## replicate(1,2*hw+1)
;       for jj=0,5 do begin
;         extract_image, profile_rect, (rect_ivar >0), xspot, sigma, $
;            f, finv, ansimage=ansimage, npoly=0,wfixed=[1,1],ymodel=ymodel
;           
;         ansimage = reform(ansimage, 2, 5, 2*hw+1)
; 
;         sigma_shift = djs_avsigclip(ansimage[1,*,*],3) / $ 
;                       (djs_avsigclip(f,1) > 1) 
;         sigma_shift = (sigma_shift < 0.5) > (-0.3) 
;
;         sigma = sigma * (1 + (sigma_shift ## replicate(1,2*hw+1)))
;
;         sigma = sigma < (min(aper))/2.
;       endfor
;
;       obj_str[q].spatial_fwhm     =  sigma[0,*] * 2.3548
;       obj_str[q].spatial_fwhm_err =  $
;           sqrt(djs_median(abs(ansimage[1,*,*]),3) / total(f,1))*2.3548


       if median_sn2 LE LOWSNR^2 then begin
          print, 'S/N is too low to recover robust profile model'
          print, 'Using other orders to choose a representative gaussian'
          other_orders = where(fwhm_here GT 0, nother)
          apodization= 0.8
          profile_cen = 0.
          cutoff_aper=aper

          if nother LT 2 then begin
               print, 'not using a fit, just taking total(aper)/12 for sigma'
               profile_fwhm = total(aper)/12. * 2 * 2.358
          endif else begin
               fwhm_coeffs = $
                  ladfit(other_orders, fwhm_here[other_orders])
               profile_fwhm = poly([q], fwhm_coeffs)
               profile_fwhm = profile_fwhm[0]

          endelse

          profile_lwhm = -0.5*profile_fwhm
          profile_rwhm =  0.5*profile_fwhm

          usesigma = profile_fwhm/2.3548

;
;           small 0.5% correction to get total area = 1?
;
          final_obj_profile = exp(-0.5*(slit_dist/usesigma)^2) $ 
              / sqrt(2.0*!Pi) / usesigma 
          profile_max =  max(final_obj_profile)

          obj_profile = final_obj_profile 

       endif else begin
           
         ;; Looks like edges of the CCD require profile changes
         do2d = (median_sn2 GT HIGHSNR)

         upper = 3 + sqrt(median_sn2)/20.
         lower = 3 + sqrt(median_sn2)/20.
         bkspace = 0.1 
         ;; Fit the object profile to the flux
         if keyword_set(do2d) then begin
             ;; High S/N case
             profile_set = bspline_iterfit(slit_dist, profile_temp, $
                                   invvar=profile_ivar, yfit=obj_profile, $
                                   maxiter=50L, upper=upper, lower=lower, $
                                   bkspace=bkspace, /groupbadpix, maxrej=2,$
                                   x2=ystart, xmin=0, xmax=nrow, npoly=3, $
                                   outmask=prof_mask, /silent) 
;                                  bkspace=0.8, /groupbadpix, maxrej=1,$
             mean_profile = bspline_valu(slit_dist, profile_set, $
                             x2=slit_dist*0 + 0.5*(profile_set.xmin + profile_set.xmax))
         endif else begin
             ;; Standard case
             profile_set = bspline_iterfit(slit_dist, profile_temp, $
                                   invvar=profile_ivar, yfit=obj_profile, $
                                   maxiter=50L, upper=upper, lower=lower, $
                                   bkspace=bkspace, /groupbadpix, maxrej=1, $
                                   outmask=prof_mask, /silent)
             mean_profile = obj_profile
         endelse

         ;; Calculate profile diagnostics
         xs = sort(slit_dist)
         profile_max = max(mean_profile[xs] * (slit_dist[xs] GE -0.5*aper[0]) * $
                                         (slit_dist[xs] LE 0.5*aper[1]) , max_place)
         left_side = where(mean_profile[xs[0:max_place]] LT profile_max/2.0)
         if left_side[0] EQ -1 then profile_lwhm = min(slit_dist) $
         else profile_lwhm = slit_dist[xs[max(left_side)]]

         right_side = where(mean_profile[xs[max_place:*]] LT profile_max/2.0)
         if right_side[0] EQ -1 then profile_rwhm = max(slit_dist) $
         else profile_rwhm = slit_dist[xs[min(right_side)+max_place]]
         profile_fwhm = profile_rwhm - profile_lwhm
         profile_cen = 0.5*(profile_rwhm + profile_lwhm)
         fwhm_here[q] =  profile_fwhm
         ;; The following returns NaN for median_sn2 < 1  
         if median_sn2 LT 1. then apodization = 0.8 else $
           apodization = x_dierfc(1./(median_sn2 / $
                                    (profile_fwhm > 0.2))^2)/2.3548 > 0.8

         if 1.7*apodization*profile_fwhm GT total(aper) then begin 
           if keyword_set(debug) then print, $
               'X_EXTECHOPT: aperture is too small, less than 2.5 FHWM, resetting'
           aper = aper * 1.7*apodization*profile_fwhm / total(aper)
         endif

         cutoff_aper = [-1.7*apodization*(profile_lwhm-profile_cen), $
                         1.7*apodization*(profile_rwhm-profile_cen)]

     endelse
     ;; Kludge
     cutoff_aper = cutoff_aper > MIN_CUT
     obj_str[q].spatial_fwhm = profile_fwhm * slit_len / 2.
     temp_profile = obj_profile / slit_proj * 2.

     ;
     ;  Check trace center here
     ;
     max_x = max(xstart)
     min_x = min(xstart)

     in = where(slit_dist - profile_cen GT 4.0*profile_lwhm AND $
                slit_dist - profile_cen LT 4.0*profile_rwhm AND $
                prof_mask GT 0, nin)

;;;;   11/11/05  Scott had some really bad code below, took care of
;;;;   that above:  min_x, max_x were hard limits to the profile
;;;;   fitting region, and not appropriate for traces near the edge of
;;;;   the slit.  We can't remember why Scott wrote the code in the
;;;;   first place, probably delusions of grandeur.
;
;     in = where(slit_dist - profile_cen GT 4.0*profile_lwhm AND $
;               xstart GT min_x +2 AND $
;                slit_dist - profile_cen LT 4.0*profile_rwhm AND $
;                xstart LT max_x - 2 AND $
;                prof_mask GT 0, nin)
  
     smooth_length = 41
     trace_corr = 0.
     if nin GT 3*smooth_length then begin
      ks = sort(ywave[in])
      k0  = smooth((temp_profile[in] * skysub[inorder[in]])[ks],  smooth_length)
      kp1 = smooth((temp_profile[in-1] * skysub[inorder[in]])[ks], $
                     smooth_length)
      km1 = smooth((temp_profile[in+1] * skysub[inorder[in]])[ks], $
                     smooth_length)

      denom = 2.0* (km1 + kp1 - 2*k0)
      denom[0:smooth_length*1.1] = 0
      denom[nin-smooth_length*1.1:*] = 0
      trace_deviation = (km1 - kp1) / (denom + (denom EQ 0)) * (denom LT 0)
     
      trace_dev_ivar = (k0 > 0)*(denom LT 0) * (abs(trace_deviation - $
                                             median(trace_deviation)) LT 1)

      trace_sn = sqrt(total(trace_dev_ivar) > 0)
      damped_corr = 0.

      if trace_sn GT 10.0 then begin

        t_nbkpt = (long(trace_sn / 40.0) < 15) > 1 
        check = where(trace_dev_ivar GT 0)
        if trace_sn GT 50 then begin
          trace_dev_set = bspline_iterfit(ywave[in[ks[check]]], $
                   trace_deviation[check], invvar = trace_dev_ivar[check], $
                   nbkpt=t_nbkpt, yfit=trace_dev_fit, upper=2, lower=2, $
                   /groupbadpix, maxrej=20, /sil,outmask=trace_dev_mask)
          trace_corr = bspline_valu(ycol, trace_dev_set)
        endif else begin
          med_trace = median(trace_deviation[check])
          trace_dev_fit = trace_deviation[check]*0.0 + med_trace
          trace_corr = ycol*0.0 + med_trace
        endelse
     
        damped_corr = trace_corr / (1. + abs((trace_corr-1.0) > 0))

        if keyword_set(chk) then begin 
          mean_trace_dev = mean(abs(trace_dev_fit))
          print, 'Possible Trace correction plotted here: ', mean_trace_dev, $
               trace_sn, t_nbkpt, format='(a,f9.5, f9.3, i6)'
          if i+i_trc EQ 0 then window, 1, title='Trace Residuals ' + obj_name  $
          else wset, 1
          plot, ywave[in[ks]], trace_deviation,ps=3, yr=[-0.5, 0.5], /xs, $
             ytitle='Pixel Deviation', xtitle='Pixel Row', $
             title='Residuals and Fit for Order'+ $
                          string(ordr_str[q].order, format='(i4)')

          oplot, ycol, trace_corr, color=clr.red
          oplot, ycol, damped_corr, color=clr.green

          if i+i_trc EQ 0 then window, 0, title='Object Profile '+ obj_name  $
          else wset, 0
          if mean_trace_dev GT 0.1 then print, 'Wow, bad trace?'
        endif
      endif else print, 'Not enough Signal to measure trace_deviation'
     endif 

     ;
     ;  Don't allow trace_corr to be larger than 1 pixel
     ;  damped_corr must fall between -1 and 1.
     ;    

     obj_str[q].fin_trc[0:nrow-1] = objcen + damped_corr

     replace_left = -1L
     replace_right = -1L
     if NOT keyword_set(final_obj_profile) then begin
     ;
     ;  Check wings of Object profile 
     ;


         final_obj_profile = obj_profile
         sort_profile = sort(obj_profile)
         low_level = obj_profile[sort_profile[n_elements(obj_profile)*0.01]]
         if low_level LT -0.05*profile_max then begin
             print, 'x_extechopt: profile is significantly negative, ', low_level
             print, '     perhaps you should issue /skipskysub option'
             print, '     or change aperture size before sky subtraction.'
         endif

            
         
         ; Fit right side :
         x = slit_dist[xs] - profile_cen
         poly_base = fpoly(cutoff_aper,3)

         fit_right = where(x GT profile_rwhm  AND x LE cutoff_aper[1] AND $
                                final_obj_profile[xs] GT 0, nright)
         replace_right = where(x GT cutoff_aper[1], nrr)
            
         if (nrr GT 0 AND nright GT 10) then begin
             right_coeff = poly_fit(x[fit_right], $
                                    alog(final_obj_profile[xs[fit_right]]), 2)
             eval_aper = total(poly_base[1,*] * right_coeff)
             eval_slope = total(poly_base[1,*] * right_coeff[1:*]*(findgen(2)+1))
             ;; JXP Kludge  10oct05
             eval_slope = eval_slope < (-1.0)
             final_obj_profile[xs[replace_right]] = $
               exp(eval_aper + $
                   (x[replace_right] - cutoff_aper[1])*eval_slope) > 0
         endif
;         print, 'Right slope', eval_aper, eval_slope
         
         fit_left = where(x LT profile_lwhm AND x GE -1.0*cutoff_aper[0] AND $
                          final_obj_profile[xs] GT 0, nleft)
         replace_left = where(x LT -1.0*cutoff_aper[0], nlr)
         if (nlr GT 0 AND nleft GT 10) then begin
             left_coeff = poly_fit(-1.0*x[fit_left], $
                                   alog(final_obj_profile[xs[fit_left]]), 2)
             eval_aper = total(poly_base[0,*] * left_coeff)
             eval_slope = total(poly_base[0,*] * left_coeff[1:*]*(findgen(2)+1))
             ;; JXP Kludge  10oct05
             eval_slope = eval_slope < (-1.0)
             final_obj_profile[xs[replace_left]] = $
               exp(eval_aper - $
                   (x[replace_left] + cutoff_aper[0])*eval_slope) > 0
         endif 
;         print, 'Left slope', eval_aper, eval_slope
     endif ; End final_object profile determination

     if keyword_set(chk) then begin 
             yrange = [-0.3, 1.5*max(final_obj_profile)]
;             xrange = slit_dist[[xs[min(where(profile_ivar[xs] GT 0)) > 0, $
;                                max(where(profile_ivar[xs] GT 0))]]]
  
           plot, slit_dist, obj_profile, /nodata, color=clr.white, $
              background=clr.black, yrange=yrange, ystyle=1, /xs, $
              xtitle='Slit Position (fraction of slit)', $
              ytitle='Object Profile Cross-section'

           oplot, slit_dist, profile_temp, ps=3, color=clr.white

           if replace_left[0] NE -1 then $
             oplot, slit_dist[replace_left], obj_profile[replace_left], $
                   ps=3, color=clr.red
           if replace_right[0] NE -1 then $
             oplot, slit_dist[replace_right], obj_profile[replace_right],$
                   ps=3, color=clr.red

           oplot, slit_dist, final_obj_profile, ps=3, color=clr.green

           oplot, -1.0*[cutoff_aper[0], cutoff_aper[0]] + profile_cen, $
                       [0.0, 1.0], thick=3, color=clr.blue
           oplot, [cutoff_aper[1], cutoff_aper[1]] + profile_cen, $
                       [0.0, 1.0], thick=3, color=clr.blue

           oplot, [-1.0*base_aper[0], -1.0*base_aper[0], $
                   base_aper[1], base_aper[1]] + profile_cen, $
                       [2.0, 0.0, 0.0, 2.0], thick=3, color=clr.burlywood

           xyouts, [0.6,0.7,0.7,0.7,0.7], [0.9,0.82,0.75,0.685,0.625],/normal, $
             [string("Order # ",ordr_str[q].order, format='(a,i3)'),$
              'Fitted Profile', 'Applied Profile', $
              'Extrapolation Limits', 'Boxcar Limits'], chars=1.5

         oplot, [profile_lwhm, profile_lwhm, profile_rwhm, profile_rwhm], $
            [profile_max/2, profile_max/2.2, profile_max/2.2, profile_max/2], $
                  thick = 3
         xyouts, [0.0], [profile_max/2.7], $
                [string(obj_str[q].spatial_fwhm, '"', format='(f6.3,a)')], $
                 alignment=0.5, charsize=1.5
 
           xr = !x.crange[1] - !x.crange[0]
           yr = !y.crange[1] - !y.crange[0]
           oplot, [0.6,0.68]*xr + !x.crange[0], $
                 [0.77,0.77]*yr + !y.crange[0], color=clr.green
           oplot, [0.6,0.68]*xr + !x.crange[0], $
                 [0.84,0.84]*yr + !y.crange[0], color=clr.red
           oplot, [0.6,0.68]*xr + !x.crange[0], $
                 [0.70,0.70]*yr + !y.crange[0], color=clr.blue
           oplot, [0.6,0.68]*xr + !x.crange[0], $
                 [0.63,0.63]*yr + !y.crange[0], color=clr.burlywood

     endif
             
      

     check = where(ivar[inorder[ys]] GT 0, nc)
     everyn = 1.0*max(histogram(ystart[check])) 

;
;     Choose velocity spacing here:
;
     if profile_rwhm LT profile_lwhm then begin
         print, 'WARNING: x_extechopt is busted, '
         print, 'Right half-maximum is less than left half-max'
         stop
     endif

     goodpix = where((wave_order GT 3) AND (ivar[inorder] GT 0) AND $
                     ((slit_dist-profile_cen) GT profile_lwhm) AND $
                     ((slit_dist-profile_cen) LT profile_rwhm), ngood)

     if ngood LT 500L then begin
         print, 'x_extechopt:  Too few good pixels in order. Skipping...'
         flg_skip = 1
         continue
     endif
     flg_skip = 0

     minlog = min(wave_order[goodpix[everyn:*]])
     maxlog = max(wave_order[(reverse(goodpix))[everyn:*]])

     if keyword_set(img_arc) then begin
;       velpix = (side EQ 1 ? 1.50d : 2.10d) * 4096./nrow 
       loglam = alog10(1.0d + velpix / 299792.458d)
       wave0  = alog10(3000.0d)
     endif else begin
       loglam = 1.0
       wave0  = 0.0d
     endelse


;;
;;    Try not to grab the very edge
;;
     npix = long((maxlog - minlog)/loglam - 1) > 2 
     loglam0 = (long((minlog - wave0)/loglam) + 2) * loglam + wave0
     logwvarr = dindgen(npix) *loglam + loglam0
     obj_str[q].npix = npix

     if keyword_set(img_arc) then obj_str[q].wave[0:npix-1] = 10.0^logwvarr $
     else obj_str[q].wave[0:npix-1] = logwvarr

     fullbkpt = [2*logwvarr[0] - reverse(logwvarr[1:bspline_ord-1]), logwvarr,$
             2*logwvarr[npix-1] - reverse(logwvarr[npix-bspline_ord:npix-2]) ]
    
     ;; Worry about errors
     if npix LT 0.6*nrow then begin
               print, 'x_extechopt:  WARNING -- You are binning too much!!'
               print, 'x_extechopt:  You should extract with ~nrow pixels '+$
                 'and then bin...unless this is a partial order.'
     endif
      

     ;;  Some helpful logging (all on one line?)
     
     print,'                                   cutoff  ---Round 1----  ---Round 2--- ' 
     print,' q Ordr Bk S/N_med Apodiz   FWHM    aper   II Chi_nu  Bad  II Chi_nu  Bad' 

     print, q, ordr_str[q].order, everyn, sqrt(median_sn2), apodization, $
       obj_str[q].spatial_fwhm, cutoff_aper, $
       format='(i2.2, i4.3, i4.2, f8.2, f7.3, f7.3, f5.1, f4.1, $)'
     
;     stop

     if obj_str[q].spatial_fwhm LT 0.4 then begin
         print, 'WARNING: spatial FWHM is less than 0.4", '
         print, '  This is either great, wrong or very undersampled '
     endif

     ;;;;;;;;;;;;
     ;; Scattered light kludge
;     if I_TRC EQ 0 then begin
;         sub = where(abs(10.^wave_order[ys]-3860.13) LT 1.)
;         img[inorder[ys[sub]]] = img[inorder[ys[sub]]] + 0.7
;     endif

     ;; Fit the flux with our optimal profile
     final_obj_profile = final_obj_profile / slit_proj * 2.
     extract_set = bspline_extract(wave_order[ys], img[inorder[ys]], $
                             ivar[inorder[ys]], slit_profile[ys], $
                             maxiter=25L, fullbkpt=fullbkpt, nord=bspline_ord, $
                             final_obj_profile[ys], yfit=yfit, $
                             /groupbadpix, $
                             maxrej=3, upper=upper,lower=lower, $
                             outmask=outmask, /relative, /silent, $
                             buff=string(' ', ' ', format='(a,39x,a)'))
;     stop

     ;; Boost the inverse variance
;     TST_IVAR = 1
;     sout = outmask[sub]
;     print, 'Mean: ', mean((img[inorder[ys[sub]]])[where(sout)])

     if keyword_set(TST_IVAR) then begin
         sub = where(abs(10.^wave_order[ys]-3860.13) LT 1.)
         spix = ys[sub]
         sout = outmask[sub]
         swv = 10.^wave_order[spix]
         chi2_spix = (img[inorder[spix]] - yfit[sub])^2 * ivar[inorder[spix]]
         chi_spix = (img[inorder[spix]] - yfit[sub]) * sqrt(ivar[inorder[spix]])
         s_reduced_chi2 = total(chi2_spix * sout)/(total(sout)- 1)
         
         sdiff= (img[inorder[spix]]-yfit[sub])[where(sout)]
         plothist, sdiff
         print, 'Standard: ', mean(sdiff), s_reduced_chi2
         
         x_splot, (img[inorder[spix]])[where(sout)], $
                  (ivar[inorder[spix]])[where(sout)], /blo, psym1=1 
         
         ;; Try again
         delvarx, outmask
         readno = 3.0
         tst_ivar = ivar
         tst_ivar[inorder[ys]] = 1.0/( (img[inorder[ys]]-yfit - 7.)>0. $
           + 7. + readno^2)

         extract_set = bspline_extract(wave_order[ys], img[inorder[ys]], $
                                       tst_ivar[inorder[ys]], slit_profile[ys], $
                                       maxiter=25L, fullbkpt=fullbkpt, $
                                       nord=bspline_ord, $
                                       final_obj_profile[ys], yfit=yfit, $
                                       /groupbadpix, $
                                       maxrej=3, upper=upper, lower=lower, $
                                       outmask=outmask, /relative, /silent, $
                                       buff=string(' ', ' ', format='(a,39x,a)'))

         chi2_spix = (img[inorder[spix]] - yfit[sub])^2 * tst_ivar[inorder[spix]]
         chi_spix = (img[inorder[spix]] - yfit[sub]) * sqrt(tst_ivar[inorder[spix]])
         s_reduced_chi2 = total(chi2_spix * sout)/(total(sout)- 1)
         
         sdiff= (img[inorder[spix]]-yfit[sub])[where(sout)]
         plothist, sdiff
         print, 'Test: ', mean(sdiff), s_reduced_chi2
         
         x_splot, (img[inorder[spix]])[where(sout)], $
                  (tst_ivar[inorder[spix]])[where(sout)], /block, psym1=1 
         
;     s_reduced_chi2 = total(chi2_spix * sout)/(total(sout)- 1)
;     chi2_spix_B = (img[inorder[spix]] - yfit[sub]+0.2)^2 * ivar[inorder[spix]]
;     s_reduced_chi2_B = total(chi2_spix_B * sout)/(total(sout)- 1)
     endif


;     gd = where(ivar[inorder[ys]] GT 0. and img[inorder[ys]] LT 15.)
;     chi_pix = (img[inorder[ys[gd]]] - yfit[gd]) * sqrt(ivar[inorder[ys[gd]]])
;     print, 'chi = ', median(chi_pix), mean(chi_pix)

     chi2_pix = (img[inorder[ys]] - yfit)^2 * ivar[inorder[ys]]
     reduced_chi2 = total(chi2_pix * outmask)/$
           (total(outmask)-n_elements(extract_set.coeff) - 1)

     tempset = create_bsplineset(extract_set.fullbkpt, extract_set.nord)
     tempset.bkmask = extract_set.bkmask
       
     ;; Create the model of the sky
     skyset = tempset
     skyset.coeff = extract_set.coeff[1,*]  ;; Coefficients for the sky
     objset = tempset
     objset.coeff = extract_set.coeff[0,*]  ;; Coefficients for the object
     isig = extract_set.icoeff[0,*]
     print, reduced_chi2
     ;; Only smooth 10-sigma rejections
     smooth_outmask = x_smoothmask(chi2_pix LT 100.0*(reduced_chi2 > 1), $
                           inorder[ys], ncol)

     if I_TRC EQ FIN_TRC then begin

         skybkpts = 0
         if skyfil[0] NE '' then begin
             skyfilename = findfile(skyfil[0]+'*', count=chkfil)
             if chkfil EQ 1 then skystrct = xmrdfits(skyfilename[0], q+1, /silent)
         endif
         
         if keyword_set(skystrct) then begin
             if N_elements(helio) EQ 0 then skybkpts = skystrct.fullbkpts $
             else begin
                 skybkpts = 10^skystrct.fullbkpt
                 airtovac, skybkpts
                 skybkpts = alog10(skybkpts) + helio
             endelse
             ;; Another round after rejecting more pixels
             bspline_indsky, wave_order[ys], img[inorder[ys]], $
                             ivar[inorder[ys]], slit_profile[ys], $
                             maxiter=25L, final_obj_profile[ys], $
                             yfit=yfit, skybkpt=skybkpts, $
                             /groupbadpix, maxrej=1, upper=upper, $
                             lower=lower, outmask=outmask, /relative,$
                             everyn=everyn, inmask=smooth_outmask, $
                             /silent, objset=objset, skyset=skyset
             
             isig = objset.icoeff[0,*]
         endif else begin
             ;; Another round after rejecting more pixels
             buff = string(' ', ' ', format='(a,55x,a)')
             print, buff, format='(a,$)'
             extract_set = bspline_extract(wave_order[ys], img[inorder[ys]], $
                                           ivar[inorder[ys]], slit_profile[ys], $
                                           maxiter=25L, final_obj_profile[ys], $
                                           yfit=yfit, /groupbadpix, maxrej=1, $
                                           lower=lower, upper=upper, alpha=alpha, $
                                           outmask=outmask, /relative, /silent, $
                                           fullbkpt=fullbkpt, nord=bspline_ord, $
                                           inmask=smooth_outmask, $
                                           covariance=covariance, buff=buff)
             
             ;; Original code
             tempset = create_bsplineset(extract_set.fullbkpt, extract_set.nord)
             tempset.bkmask = extract_set.bkmask
             
             ;; Create the model of the sky
             skyset = tempset
             skyset.coeff = extract_set.coeff[1,*]  ;; Coefficients for the sky
             objset = tempset
             objset.coeff = extract_set.coeff[0,*]  ;; Coefficients for the object
             isig = extract_set.icoeff[0,*]
         endelse
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
         ;; Test the region
;         sub = where(abs(10.^wave_order[ys]-3860.73) LT 0.12)
;         sout = outmask[sub]
;         sdiff= (img[inorder[ys[sub]]]-yfit[sub])[where(sout)]
;         plothist, sdiff
;         print, 'Standard: ', mean(sdiff)
;       
;         stop
     endif


       objset_err = objset
       objset_err.coeff = 1./(isig + (isig EQ 0)) * (isig GT 0)
       sky_spectra  = bspline_valu(wave_order, skyset)
       if zero_wave[0] NE -1 then sky_spectra[zero_wave] = 0.
       model_sky[inorder]    =  sky_spectra * slit_profile

       ;; Update sky_guess for profile determination
       sky_guess = skysub[inorder] - img[inorder] + model_sky[inorder]

       ;; Calculate the final flux model 
       obj_spectra  = bspline_valu(wave_order, objset)
       if zero_wave[0] NE -1 then obj_spectra[zero_wave] = 0.
       model_obj[inorder]    =  obj_spectra*final_obj_profile
       model_prof[inorder]    =  final_obj_profile

;      if ordr_shift[q].order EQ 85 then stop
       ;; Mask
;       if I_TRC EQ FIN_TRC then stop


   ENDFOR  ; End trace correction and profile fitting

;   stop

   ;; Here ENDS boxcar
   if keyword_set( BOXONLY ) then continue
   if FLG_SKIP EQ 1 then continue
   
   ;;
   ;;     Calculate the 1d spectra
   ;;
   action=0
   loweraction=0
   upperaction=0
   obj_str[q].fx[0:npix-1] = bspline_valu(logwvarr, objset, action=action, $
                                          lower=loweraction, upper=upperaction)
   
   print, 'Calculating Variance', format='(a,$)'
   total_var = bspline_calcvar(action, loweraction, upperaction, alpha)
   print, '..Done'
   
   obj_str[q].var[0:npix-1] = total_var * (total_var GT 0)
   
   print, npix, obj_str[q].wave[[0,npix-1]], long(total(outmask EQ 1)), $
     format='(i10, f15.5, f15.5, i10)'
   if keyword_set( OCHK ) then begin
       x_splot, obj_str[q].wave[0:npix-1], obj_str[q].fx[0:npix-1], $
         ytwo=obj_str[q].sig[0:npix-1], /block
   endif

   ;; JXP  1/19/06  (smoothed)
;   extract_mask[inorder[ys]] = outmask
   extract_mask[inorder[ys]] = x_smoothmask(outmask, inorder[ys], ncol)

   chi_image[inorder[ys]] = (img[inorder[ys]] - yfit) * $
     sqrt(ivar[inorder[ys]] * outmask)
   
   if keyword_set( DEBUG ) then begin
       if keyword_set(img_arc) then $
         xatv, (img-(model_obj+model_sky)), wvimg=10^img_arc, $
         /block, min=-30., max=30. else $
             xatv, (img-(model_obj+model_sky)), /block, min=-30., max=30. 
   endif

;    gd = where(ivar[inorder[ys]] GT 0. and img[inorder[ys]] LT 15.)
;    chi_pix = (img[inorder[ys[gd]]] - yfit[gd]) * sqrt(ivar[inorder[ys[gd]]])
;    print, 'chi = ', median(chi_pix), mean(chi_pix)
;    model_obj[inorder]    =  obj_spectra*final_obj_profile
;    stop

;   stop

  endfor                          ; End of order loop

   ;; Display residuals
   if keyword_set( RESCHK ) then begin
       xatv, (img-(model_obj+model_sky)), /block, min=-30., max=30.
;       xatv, (skysub-(model_obj+model_sky))*sqrt(ivar), /block, min=-5., max=5.
;       xatv, (skysub-(model_obj+model_sky))*sqrt(ivar)*(extract_mask EQ 0), $
;         /block, min=-5., max=5.
   endif
;   stop



;  DONE
  print, 'x_extechopt: All done! '
  return
end
