;+ 
; NAME:
; x_slitflat
;     Version 1.1
;
; PURPOSE:
;    Stores slit profile and gradient along each order from twilight flats.
;    This routine is critical for performing ideal sky subtraction,
;    especially given the short slit length.  The following steps are
;    performed in x_slitflat_work:
;    
;    1.  Fit and subtract the scattered light in the twilight flat
;    2.  Extract a boxcar profile down the center of the flat
;    3.  Loop on orders
;    4.  Calculate the slit angle at each pixel in the order
;    5.  Calculate the Jacobian   (DEPRECATED)
;    6.  Fit a bspline to the profile
;    7.  Run diagnsotics on the fit
;    8.  Save the good ones to profile0 and profile1 tags
;
; CALLING SEQUENCE:
;   
;  x_slitflat, x, setup, side, [/chk, /clobber]
;
; INPUTS:
;   x     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;   Fills in the profile0 and profile1 tags in the Order structure
;
; OPTIONAL KEYWORDS:
;  /CHK  - Show the profiles and fits order by order
;  /CLOBBER - Clobber an previous work
;  TFLAT_FIL - Name for TFLAT file
;  RESIDUAL_FIL  - Output name for Jacobian matrix
;  PROFILE_FIL - Output name for profile fits
;  DETILT    - Remove a linear tilt from the Cross-section fit
;
;  NXBKT  -  Number of x breakpoints for scattered light fit 
;  NYBKT  -  Number of x breakpoints for scattered light fit 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_slitflat, mike, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_ordermask
;  x_modelslit
;  x_qw
;  x_slitflat_work
;
; REVISION HISTORY:
;   ??-??-2003 Written by SB
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_slitflat, tflat, tflativar, ordr_str, $
                  maskimage=maskimage, tflat_file=tflat_file, gapfit=gapfit, $
                  chk=chk, residual=residual, prof_image=prof_image, $
                  qa_str=qa_str, detilt=detilt, NOSCATT=noscatt, $
                  tflat_spec_fil=tflat_spec_fil, SLIT_SAMP=slit_samp

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'fil = x_slitflat(tflat, tflativar, ordr_str, /clobber, /chk, /DETILT' 
      print, '    RESIDUAL_FIL=, TFLAT_FIL=, PROFILE_FIL=  [v1.1]'
      return,-1
  endif

  if NOT keyword_set(slit_samp) then slit_samp = 0.06
  
  if keyword_set(tflat_file) then begin        
      tflat     = xmrdfits(tflat_file,0)
      tflativar = xmrdfits(tflat_file,1)
  endif

;    Colors
  clr = getcolor(/load)
  
  sz        = size(tflat, /dimen)
  nord = n_elements(ordr_str)
  tflat_spec = fltarr(sz[1], nord)
  
;
;    Mapping out orders and gaps, this is kind of slow
;
  print, 'Calling x_ordermask...', format='(a,$)'
  maskimage = x_ordermask(sz[0], sz[1], ordr_str, trim=0.7)
  print, 'Done.'
  
  ncol= (size(maskimage))[1]
  nrow= n_elements(ordr_str[0].lhedg)
  
;  if NOT keyword_set(nxbkpt) then nxbkpt=5
;  if NOT keyword_set(nybkpt) then nybkpt=10
  
  t0 = systime(1)
;  x = findgen(sz[0])
;  y = findgen(sz[1])
  
  if not keyword_set( NOSCATT ) then begin
      print, 'x_slitflat: Trying to fit scattered light background'
;      print, 'x_slitflat:  Fitting Interorder Light...', nxbkpt, nybkpt, format='(a,i4,i4,$)'
      x_modelslit, tflat, tflativar, ordr_str, scat_model=gapfit
;      model_slit, tflat, tflativar, ordr_str, scat_model=gapfit
      tflat_sub = tflat  - transpose(x_medianrow(transpose(gapfit), 31))
  endif else tflat_sub = tflat
          
  residual = tflat * 0.0
  prof_image = tflat * 0.0
  
  
  profile_x        = findgen(251)/100. - 1.25
  if NOT keyword_set(focus) then focus = 11.5
;   oprof = 0.5 - errorf(focus*(abs(profile_x)-1))/2.0
  oprof = 0.0
  
  ordr_str.profile0 = 0
  ordr_str.profile1 = 0
  qa_str = 0
  
  if keyword_set(chk) then $
    print, ' i Ord   N_ft   N_1d   N_2d  linear_coeffs  ----chi_1d----   ----chi_2d----'

 
   ;; Loop on the Orders
   for iord=min(ordr_str.order), max(ordr_str.order),1 do begin
      q = where(ordr_str.order EQ iord)
      q = q[0]
      if q EQ -1 then continue

      if ordr_str[q].flg_anly EQ 0 then begin
         print, 'x_slitflat: Skipping order', ordr_str[q].order
         continue
      endif

      inorder = where(maskimage EQ iord  AND tflativar GT 0, nin)
;      inorder = where((maskimage EQ iord  OR $
;                       maskimage EQ -1*iord OR $
;                       maskimage EQ -1*iord +1) AND $
;                      tflativar GT 0, nin)
      print, q, iord, format='(i2, i4, $)'


       
     
      if nin LE nrow then begin
          print, '...not enough pixels in order'
          continue
      endif 

      ystart = 1.0d*(inorder / ncol)
      oo = (ordr_str[q].lhedg + ordr_str[q].rhedg)/2.0
      ordrcen =  oo[ystart] 
      xstart = 1.0d*(inorder mod ncol)  
      ycol = dindgen(nrow)
      y = ((2*ystart - nrow) / nrow)

      ;; Calculate wavelengths (and the slope at each x,y pair)
      ywave = x_qckwav(xstart - ordrcen, ystart, ordr_str[q].arc_m, $
                        arc_slope=arc_slope)

      slit_frac = frac_order(ordr_str[q], xstart, ywave)
      slit_proj = slitproj(ystart, ywave, arc_slope, ordr_str[q])

      ;;
      ;;  now estimate flux_tflat(lambda)
      ;;     
      N = n_elements(inorder)
    
      xsort = sort(slit_frac)
      tflat_norm = tflat_sub[inorder]  
;
;    good will NOT work on dual slit... we need to find a way to get
;    pixel
;
      good = where(abs(slit_frac) LT 0.7 AND $
                   xstart GT 4 AND xstart LT ncol-4 AND $
                   tflativar[inorder] GT 0,ngood)

      maxy = max(ywave[good])
      miny = min(ywave[good]) 
      everyn = long(1.2*ngood/(maxy - miny))
      ys = sort(ywave[good])
      bset = bspline_iterfit(ywave[good[ys]], tflat_norm[good[ys]], $
                      everyn=everyn, /groupbadpix, upper=10, lower=10, $
                      invvar=tflativar[inorder[good[ys]]], /silent)
      tflat_fit = bspline_valu(((ywave < maxy) > miny), bset)
      tflat_spec[*,q] = bspline_valu(ycol, bset)

;
;  Now it's time for the slit profile
;

      profile = tflat_norm / (tflat_fit + (tflat_fit EQ 0)) * (tflat_fit GT 20)

      profile_sub = profile 

      profile_ivar = ((tflativar[inorder] * tflat_fit^2) < 1.0e5)  * $
        (profile GT -0.1 AND profile LT 1.5 AND tflat_fit GT 20)


      profile_med = median(profile[xsort],15)
      outlier = (profile[xsort] - profile_med)^2*profile_ivar[xsort] LT 100
      totalivar = total(profile_ivar[xsort]*outlier) / total(outlier)
      slit_samp_use = slit_samp 
     
      bkpts0 = 0 
      profile_set = bspline_iterfit(slit_frac[xsort], profile_sub[xsort], $
                                    invvar=profile_ivar[xsort]*outlier, $
                                    upper=2, lower=2, $
                                    /groupbadpix, maxrej=10, bkpt=bkpts0, $
                                    bkspace=slit_samp_use, yfit=profile_fit, $
                                    outmask=profile_mask, /silent) 

      if size(profile_set,/tname) EQ 'INT' then begin
          print, 'Could not fit 1d bspline'
          continue
      endif

      if keyword_set(chk) then plot, ywave, tflat_fit, ps=3, /xs

      ;; JXP: Above fit looks a little low? 
      ;; Let's try just a slope fit, with different bkpts...
      ;; All fits so far look like we could independently fit for the slope.

      profile_res = (profile_sub[xsort] - profile_fit) / (y[xsort]+(y[xsort] EQ 0))
      profile_res_ivar = (profile_ivar[xsort]) * y[xsort]^2
      ;; This next stuff is to set the breakpoints and 
      ;; deal with significant order overlap
      inner = where(abs(bkpts0) LT 0.8, n_inner)
      outer = where(abs(bkpts0) GE 0.8, n_outer)
      if n_outer NE 0 then begin
          bkpts1 = bkpts0[where(abs(bkpts0) GE 0.8)] 
          if n_inner GE 4 then $
            bkpts1 = [bkpts1, $
                      bkpts0[inner[2*lindgen(n_inner/2-1)+1]]]
      endif else bkpts1 = bkpts0[inner[2*lindgen(n_inner/2-1)+1]]
      bkpts1 = bkpts1[sort(bkpts1)]

      profile_setslope = bspline_iterfit(slit_frac[xsort], profile_res, $
                                 invvar=profile_res_ivar*outlier, $
                                 bkpt=bkpts1, yfit=profile_fit_slope, $
                                 upper=2, lower=2, $
                                 /groupbadpix, maxrej=10, $
                                 outmask=profile_mask2d, /silent) 
      
      if size(profile_set2d,/tname) EQ 'INT' then begin
          print, 'Could not fit 2d bspline'
          continue
      endif

;      profile_fit = profile_fit + orig_profile[xsort]
      profile_fit2d = profile_fit + profile_fit_slope*y[xsort]

      scatter1d = total(abs((profile[xsort] - profile_fit) *profile_mask) +$
                        (1-profile_mask))/ (total(profile_mask)+1)
      scatter2d = total(abs((profile[xsort] - profile_fit2d)*profile_mask2d)+$
                        (1-profile_mask2d))/ (total(profile_mask2d)+1)
      
      print, scatter1d, scatter2d, totalivar, slit_samp_use, n_elements(bkpts1)
     
      if (totalivar LT 1000.0) then begin
          if (totalivar LT 200.0) then begin
              print, 'x_slitflat: Not enough signal to fit anything'
              print, 'x_slitflat: Skipping for now'
          endif else begin
              print, 'x_slitflat: Not enough signal to fit 2d'

              ordr_str[q].profile0 = $
                (oprof + bspline_valu(profile_x, profile_set)) * $
                (profile_x LT max(slit_frac[xsort]*profile_mask) AND $
                 profile_x GT min(slit_frac[xsort]*profile_mask))
              prof_image[inorder[xsort]] = profile_fit * tflat_fit[xsort]
              residual[inorder[xsort]] = (profile[xsort] - profile_fit)* $
                sqrt(profile_ivar[xsort] * profile_mask) 
          endelse
      endif else begin
          
;
;    This is the linear change with row  (over nrow/2)
;
          ordr_str[q].profile0  = $
            (oprof + bspline_valu(profile_x, profile_set)) * $
            (profile_x LT max(slit_frac[xsort]*profile_mask2d) AND $
             profile_x GT min(slit_frac[xsort]*profile_mask2d))
          ordr_str[q].profile1  = bspline_valu(profile_x, profile_setslope) * $
            (profile_x LT max(slit_frac[xsort]*profile_mask2d) AND $
             profile_x GT min(slit_frac[xsort]*profile_mask2d))
          
          prof_image[inorder[xsort]] = profile_fit2d * tflat_fit[xsort]
          residual[inorder[xsort]] = (profile[xsort] - profile_fit2d)* $
            sqrt(profile_ivar[xsort] * profile_mask2d) 
      endelse

      djs_iterstat, (profile[xsort] - profile_fit2d)* $
        sqrt(profile_ivar[xsort]), sigma=sig2d, sigrej=5.
      djs_iterstat, (profile[xsort] - profile_fit2d)* $
        sqrt(profile_ivar[xsort]*profile_mask2d),sigma=sig2d_clean, sigrej=5.
      djs_iterstat, (profile[xsort] - profile_fit)* $
        sqrt(profile_ivar[xsort]), sigma=sig, sigrej=5.
      djs_iterstat, (profile[xsort] - profile_fit)* $
        sqrt(profile_ivar[xsort]*profile_mask), sigma=sig_clean, sigrej=5.
    
      x = slit_frac[xsort] 
      ab = func_fit(x, profile_fit, 2, $
           invvar=profile_ivar*(abs(x) LT 0.7), function_name='poly')

      print, total([[outlier],[profile_mask],[profile_mask2d]],1), $
        [ab, sig,sig_clean, sig2d, sig2d_clean], format='(3(i7),6(f8.4))'
      
      temp_qa = x_mkslitqa(slit_frac[xsort], profile[xsort], profile_mask, $
                     ordr_str[q])
      temp_qa.ab = ab
      temp_qa.npix = n_elements(profile_mask)
      temp_qa.nrej = total(profile_mask EQ 0)
      temp_qa.chi2 = sig_clean^2

      qa_str = struct_append(qa_str, temp_qa)

      if keyword_set(chk) then begin 
          plot, slit_frac, profile, ps=3, yr=[0.9, 1.1], color=clr.white, $
            background=clr.black, xr=[-1.3,1.3], /xs
          oplot, slit_frac[xsort], profile_fit2d, ps=3, color=clr.red
          oplot, slit_frac[xsort], profile_fit, ps=3, color=clr.blue
          oplot, profile_x, ordr_str[q].profile0, color=clr.green
          oplot, [-1.3,1.3], poly([-1.3,1.3],ab), color=clr.orange, thick=3
          
          oplot, slit_frac, profile+1, ps=3
          oplot, slit_frac[xsort], profile_fit2d+1, ps=3, color=clr.red
          oplot, slit_frac[xsort], profile_fit+1, ps=3, color=clr.blue
          oplot, profile_x, ordr_str[q].profile0+1, color=clr.green
      endif

      ;; remove slope

      if keyword_set(detilt) then begin
        ss = poly(profile_x, ab)
        ordr_str[q].profile0 = ordr_str[q].profile0 / (ss + (ss EQ 0)) * $
                    (ss GT 0.5) 
        ordr_str[q].profile1 = ordr_str[q].profile1 / (ss + (ss EQ 0)) * $
                    (ss GT 0.5) 
      endif

  endfor

  if keyword_set(tflat_spec_fil) then begin
    print, 'Storing TFLAT spectra in ', tflat_spec_fil
    mwrfits, tflat_spec, tflat_spec_fil, /create
    if keyword_set(GAPFIT) then mwrfits, gapfit, tflat_spec_fil
  endif


;
;    Fix any bad orders
;
  ordermask  = total(ordr_str.profile0,1) GT 0
  badorders = where(ordermask EQ 0, nbad)
  
  if nbad GT 0 then begin
      print, 'x_slitflat: Fixing profiles with extrapolation: ', $
        strcompress(strjoin(string(badorders)))
      y0 = transpose(ordr_str.profile0)
      y1 = transpose(ordr_str.profile1)
      x = lindgen(nord) # replicate(1,251)
      ymask = (y0 NE 0)
      xy2traceset, x, y0, tset0, ncoeff=2, yfit=yfit0, invvar=ymask
      ordr_str[badorders].profile0 = transpose(yfit0[badorders,*])
      qa_str[badorders].mid = ordr_str[badorders].profile0
      qa_str[badorders].edgel = ordr_str[badorders].profile0
      qa_str[badorders].edger = ordr_str[badorders].profile0
   
  endif

  return, ordr_str
end

