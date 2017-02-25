;+ 
; NAME:
;  mike_slitflat4x2
;     Version 1.0
;
; PURPOSE:
;    Stores slit profile and gradient along each order from twilight flats.
;        Takes TFLAT as the default
;
; CALLING SEQUENCE:
;   
;  mike_mktflat, mike, setup, side
;
; INPUTS:
;   mike     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;   A new ordr_struct stored in "Flats/OStr_R_01wP.fits'
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_mktflat, mike, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Feb-2003 Written by SB
;   18-Apr-2003 Revised by JXP
;   30-Apr-2003 Polished by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mike_slitflat_work, tflat, tflativar, ordr_str, $
           maskimage=maskimage, tflat_file=tflat_file, gapfit=gapfit, $
           chk=chk, jacobian=jacobian, residual=residual

     if NOT keyword_set(slit_samp) then slit_samp = 0.04

     if keyword_set(tflat_file) then begin        
       tflat     = xmrdfits(tflat_file,0)
       tflativar = xmrdfits(tflat_file,1)
     endif

;    Colors
     clr = getcolor(/load)

     sz        = size(tflat, /dimen)
     n_orders = n_elements(ordr_str)

;
;    Mapping out orders and gaps, this is kind of slow
;
     print, 'Calling mike_ordermask...', format='(a,$)'
     wait,1. 
     maskimage = mike_ordermask(sz[0], sz[1], ordr_str, trim=1.5)
     print, 'Done.'

;
;    Fit scattered light in order gaps
;
     gapfit = mike_fitgap(tflat, tflativar, maskimage)

     tflat_sub = tflat - gapfit
;     trans_tflat = transpose(tflat - gapfit)
;     trans_ivar = transpose(tflativar)

;
;    What's up with this jacobian
;
;       |                  |
;       |  dx/ds    dy/ds  |
;       |                  |
;       |  dx/dl    dy/dl  |
;       |                  |
;
;

   chebypoly = n_elements(ordr_str[0].lhc)
   ncol= (size(maskimage))[1]
   nrow= n_elements(ordr_str[0].lhedg)
   nord = n_elements(ordr_str)
   jacobian = tflat * 0.0
   residual = tflat * 0.0
;   new_struct = struct_addtags(ordr_str, $
;         replicate({ profile0: fltarr(251), profile1: fltarr(251)},nord))

   profile_x        = findgen(251)/100. - 1.25

   for q=0,nord-1 do begin

      inorder = where(maskimage EQ ordr_str[q].order AND tflativar GT 0, nin)
      if nin GT nrow then begin

        ystart = 1.0d*(inorder / ncol)
        ordrcen =  (ordr_str[q].lhedg[ystart] + ordr_str[q].rhedg[ystart])/2.0
        xstart = 1.0d*(inorder mod ncol)  - ordrcen
        ycol = dindgen(sz[1])
        slope_spline = spl_init(ycol, double(ordr_str[q].arc_m))

        ywave = mike_qw(xstart, ystart, ordr_str[q].arc_m, arc_slope=arc_slope)


        y = ((2*ystart - nrow) / nrow)
        slit_length = ordr_str[q].rhedg[ystart] - ordr_str[q].lhedg[ystart]
        slit_length[*] = ordr_str[q].rhedg[nrow/2] - ordr_str[q].lhedg[nrow/2]
        slit_frac   = 2.0* xstart / slit_length
        test_deriv = fchebyshev_deriv1(y,chebypoly, /ideriv) 
        test_coeff = ((ordr_str[q].lhc + ordr_str[q].rhc) + $
                     (slit_frac * (ordr_str[q].rhc-ordr_str[q].lhc)))/ nrow
        trace_deriv = test_deriv # test_coeff

        ;; Jacobian
        dydl = 1.0/sqrt(1.0 + trace_deriv^2)
        dxdl = trace_deriv * dydl

   
        dxds = 1.0/sqrt(1.0 + arc_slope^2)
        dyds = arc_slope* dxds


        jacobian_ord = dxds * dydl - dxdl * dyds
        jacobian[inorder] = jacobian_ord


     ;
     ;  now estimate flux_tflat(lambda)
     ;     
     N = n_elements(inorder)
    
     xsort = sort(slit_frac)
;     tflat_norm = tflat_sub[inorder] / jacobian_ord
     tflat_norm = tflat_sub[inorder] 

;
;    good will NOT work on dual slit... we need to find a way to get
;    pixel
;
     xedge = xstart + ordrcen
     good = where(abs(slit_frac) LT 0.7 AND $
                        xedge GT 4 AND xedge LT ncol-4,ngood)

     everyn = long(0.6*ngood/(max(ywave) - min(ywave)))
     everyn = 6L  ; hard wiring for 4x2 data -- col binning is issue 
                  ; (value was 3, too small)   2x2 , 12 is a good value.
                  ; everyn is the number of pixels per free parameter.
                  ; perfect data, no rejection, you might be able to go 
                  ; down to nyquist sampling. but when rejection is needed
                  ; sampling needs to be larger.

     bset = bspline_iterfit(ywave[good], tflat_norm[good], everyn=everyn, $
                        /groupbadpix, maxrej=3, invvar=tflativar[good], /silent)
     tflat_fit = bspline_valu(ywave, bset)

     ;
     ;  Now it's time for the slit profile
     ;

     profile = tflat_norm / (tflat_fit + (tflat_fit EQ 0)) * (tflat_fit GT 0)
;     profile_ivar = tflativar * jacobian_ord^2 * tflat_fit^2 * $
;                    (profile GT -0.1 AND profile LT 1.5 AND tflat_fit GT 0)
     profile_ivar = tflativar * tflat_fit^2 * $
                    (profile GT -0.1 AND profile LT 1.5 AND tflat_fit GT 0)
 
     print, 'Working on Order...', q, ordr_str[q].order
     wait, 1.0

     profile_med = median(profile[xsort],15)
     outlier = (profile[xsort] - profile_med)^2*profile_ivar[xsort] LT 100
     profile_set = bspline_iterfit(slit_frac[xsort], profile[xsort], $
            invvar=profile_ivar[xsort]*outlier, $
            bkspace=slit_samp, yfit=profile_fit, $
            outmask=profile_mask, /silent) 

     ;; JXP: Above fit looks a little low? 

     profile_set2d = bspline_iterfit(slit_frac[xsort], profile[xsort], $
            invvar=profile_ivar[xsort]*outlier, x2=y[xsort], npoly=2, $
            bkspace=slit_samp, yfit=profile_fit2d, $
            outmask=profile_mask2d, /silent) 

     ordr_str[q].profile0  = bspline_valu(profile_x, profile_set) * $
                         (profile_x LT max(slit_frac[xsort]*profile_mask) AND $
                          profile_x GT min(slit_frac[xsort]*profile_mask))

     profile_setslope = profile_set
;
;    This is the linear change with row  (over nrow/2)
;
     profile_setslope.coeff = profile_set2d.coeff[1,*]
     ordr_str[q].profile1  = bspline_valu(profile_x, profile_setslope) * $
                         (profile_x LT max(slit_frac[xsort]*profile_mask) AND $
                          profile_x GT min(slit_frac[xsort]*profile_mask))
     
   
     residual[inorder[xsort]] = (profile[xsort] - profile_fit2d)* $
          sqrt(profile_ivar[xsort] * profile_mask2d) 

       djs_iterstat, (profile[xsort] - profile_fit)* $
          sqrt(profile_ivar[xsort]), sigma=sig, sigrej=5.
       djs_iterstat, (profile[xsort] - profile_fit)* $
          sqrt(profile_ivar[xsort]*profile_mask), sigma=sig_clean, sigrej=5.
       djs_iterstat, (profile[xsort] - profile_fit2d)* $
          sqrt(profile_ivar[xsort]), sigma=sig2d, sigrej=5.
       djs_iterstat, (profile[xsort] - profile_fit2d)* $
          sqrt(profile_ivar[xsort]*profile_mask2d), sigma=sig2d_clean, sigrej=5.

     if keyword_set(chk) then begin 
       print, [1.0*q, ladfit(slit_frac[xsort[good]], profile_fit[good]),$
                sig,sig_clean, sig2d, sig2d_clean]
       wait, 1.0
       plot, slit_frac, profile, ps=3, yr=[0.9, 1.1], color=clr.white, $
         background=clr.black
       oplot, slit_frac[xsort], profile_fit2d, ps=3, color=clr.red
       oplot, slit_frac[xsort], profile_fit, ps=3, color=clr.blue
       oplot, profile_x, ordr_str[q].profile0, color=clr.green

       oplot, slit_frac, profile+1, ps=3
       oplot, slit_frac[xsort], profile_fit2d+1, ps=3, color=clr.red
       oplot, slit_frac[xsort], profile_fit+1, ps=3, color=clr.blue
       oplot, profile_x, ordr_str[q].profile0+1, color=clr.green
      endif
    endif
     print, '           Done with that one'

    endfor

    return, ordr_str
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_slitflat, mike, setup, side, chk=chk, $
            tflat_fil=tflat_fil, ordr_fil=ordr_fil, profile_fil=profile_fil, $
            jacobian_fil=jacobian_fil, mflat_fil=mflat_fil, CLOBBER=clobber

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_slitflat, mike, setup, [side]'
      return
  endif

   if not keyword_set( SIDE ) then side = [2L]
   if not keyword_set( chk  ) then chk=0
   if not keyword_set( SLITARC  ) then slitarc = 0L

; Setup
   if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 

;  Loop on side
   for ii=0L,n_elements(side)-1 do begin
       qq = side[ii]
       ;; SIDE
       if qq EQ 1 then begin
           print, 'mike_slitflat: Profiling BLUE trace flat'
           nm = 'B'
           stop, 'Blue side is non-functional for now'
       endif else begin
           print, 'mike_slitflat: Profiling RED trace flat'
           nm = 'R' 
       endelse

       ;; Order structure
       ordr_str   = mike_getfil('ordr_str', setup, SIDE=qq)

       ;; Process an arc if necessary!
       if ordr_str[0].arc_m[0] EQ 0. then begin
           ;; Find the arc
           print, 'mike_slitflat: Processing an Arc.  Wish me luck!'
           wait, 1.
           arc = where(mike.type EQ 'ARC' AND mike.flg_anly NE 0 AND $
                       mike.setup EQ setup AND mike.side EQ qq, narc)
           if narc EQ 0 then stop else arc = arc[slitarc]
           ;; Process
           resolve_routine, 'mike_allarc', /no_recompile, /either
           rslt = mike_allarc_sngl(mike[arc].rootpth+mike[arc].img_root, $
                                   setup, qq, CLOBBER=clobber)
           if rslt EQ -1 then stop
           ;; Reread order structure and check
           ordr_str   = mike_getfil('ordr_str', setup, SIDE=qq)
           if ordr_str[0].arc_m[0] EQ 0. then stop
       endif
           
       ;; TFLAT
       tflat     = mike_getfil('tflat_fil',setup, SIDE=qq)
       tflativar = mike_getfil('tflat_fil',setup, SIDE=qq, INDX=1)

       ;; MFLAT
       mflat = mike_getfil('mflat_fil', setup, SIDE=qq)

       ;;
       print, 'mike_slitflat: Applying MFlat...'
       wait,1. 
       tflat = tflat / (mflat + (mflat EQ 0))
       tflativar = tflativar * mflat^2 * (mflat GT 0.5)


       if NOT keyword_set(profile_fil) then $
         profile_fil = mike_getfil('ordr_str', setup, SIDE=qq, /name)

       if NOT keyword_set(proimg_fil) then $
         proimg_fil = 'Flats/Profile_'+nm+'_'+c_s+'.fits'

       if NOT keyword_set(jacobian_fil) then $
         jacobian_fil = 'Flats/Jacobian_'+nm+'_'+c_s+'.fits'

       profile_struct = mike_slitflat_work(tflat, tflativar, ordr_str, $
                                           jacobian=jacobian, chk=chk, $
                                           residual=residual)

       ;; Write to fits
       print, 'mike_slitflat: Writing order w/Profile structure ', profile_fil
       wait,1. 
       mwrfits, profile_struct, profile_fil, /create
       
       ;; Write to fits
       print, 'mike_slitflat: Writing Profile Image structure ', proimg_fil
       wait,1. 
       mwrfits, profile_img, proimg_fil, /create
       
       if keyword_set(jacobian) then begin
         print, 'mike_slitflat: Writing Jacobian Image', jacobian_fil
         wait,1. 
         mwrfits, jacobian, jacobian_fil, /create
         mwrfits, residual, jacobian_fil
       endif

   endfor

   ;; All done
   print, 'mike_slitflat: All done'
   return

end


