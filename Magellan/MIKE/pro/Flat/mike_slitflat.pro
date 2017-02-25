;+ 
; NAME:
; mike_slitflat
;     Version 1.1
;
; PURPOSE:
;    Stores slit profile and gradient along each order from twilight flats.
;    This routine is critical for performing ideal sky subtraction,
;    especially given the short slit length.  The following steps are
;    performed in mike_slitflat_work:
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
;  mike_slitflat, mike, setup, side, [/chk, /clobber]
;
; INPUTS:
;   mike     -  MIKE structure
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
;  NODETILT    - Do not remove a linear tilt from the Cross-section fit
;
;  NXBKT  -  Number of x breakpoints for scattered light fit 
;  NYBKT  -  Number of x breakpoints for scattered light fit 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_slitflat, mike, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_ordermask
;  mike_fitgap
;  mike_qw
;  mike_slitflat_work
;
; REVISION HISTORY:
;   ??-??-2003 Written by SB
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;function mike_slitflat_work, tflat, tflativar, ordr_str, $
;         maskimage=maskimage, tflat_file=tflat_file, gapfit=gapfit, $
;         chk=chk, residual=residual, prof_image=prof_image, $
;;         nxbkpt=nxbkpt, nybkpt=nybkpt, qa_str=qa_str, nodetilt=nodetilt, $
;         tflat_spec_fil=tflat_spec_fil

;     if NOT keyword_set(slit_samp) then slit_samp = 0.09
;
;     if keyword_set(tflat_file) then begin        
;       tflat     = xmrdfits(tflat_file,0)
;       tflativar = xmrdfits(tflat_file,1)
;     endif

;    Colors
;     clr = getcolor(/load)
;
;     sz        = size(tflat, /dimen)
;     nord = n_elements(ordr_str)
;     tflat_spec = fltarr(sz[1], nord)

;
;    Mapping out orders and gaps, this is kind of slow
;
;     print, 'Calling mike_ordermask...', format='(a,$)'
;     maskimage = mike_ordermask(sz[0], sz[1], ordr_str, trim=0.7)
;     print, 'Done.'

;     ncol= (size(maskimage))[1]
;     nrow= n_elements(ordr_str[0].lhedg)

;     if NOT keyword_set(nxbkpt) then nxbkpt=5
;     if NOT keyword_set(nybkpt) then nybkpt=10
;     print, 'Fitting Interorder Light...', nxbkpt, nybkpt, format='(a,i4,i4,$)'

;     t0 = systime(1)
;     gapfit = tflat* 0.
;     x = findgen(sz[0])
;     y = findgen(sz[1])
;     min_gap = min(maskimage)
;     max_gap = max(maskimage[where(maskimage LT 0)])
;     scat_col = 0
;     scat_order = 0
;     scat_row = 0
;     scat_val = 0
     

;     print, 'Trying to fit scattered light background'
;     
;     model_slit, tflat, tflativar, ordr_str, scat_model=gapfit
;     tflat_sub = tflat  - gapfit
;     
;     residual = tflat * 0.0
;     prof_image = tflat * 0.0
;     
;     
;  return, ordr_str
;end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_slitflat, mike, setup, side, chk=chk, $
            tflat_fil=tflat_fil, profile_fil=iprofile_fil, $
            residual_fil=iresidual_fil, CLOBBER=clobber, NODETILT=nodetilt, $
            tflat_spec_fil=itflt_spec_fil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_slitflat, mike, setup, [side], /clobber, /chk, /NODETILT' 
      print, '    RESIDUAL_FIL=, TFLAT_FIL=, PROFILE_FIL=  [v1.1]'
      return
  endif

   if not keyword_set( SIDE ) then side = [1L,2L]
   if not keyword_set( chk  ) then chk=0

; Setup
   if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 

;  Loop on side
   for ii=0L,n_elements(side)-1 do begin
       qq = side[ii]
       ;; SIDE
       if qq EQ 1 then begin
           print, 'mike_slitflat: Profiling BLUE trace flat'
           nm = 'B'
       endif else begin
           print, 'mike_slitflat: Profiling RED trace flat'
           nm = 'R' 
       endelse

       if NOT keyword_set(iprofile_fil) then $
         profile_fil = mike_getfil('ordr_str', setup, SIDE=qq, /name) $
       else profile_fil = iprofile_fil

       if NOT keyword_set(iproimg_fil) then $
         proimg_fil = 'Flats/Profile_'+nm+'_'+c_s+'.fits' $ 
       else proimg_fil = iproimg_fil
       if x_chkfil(proimg_fil+'*',/silent) NE 0 $
         and not keyword_set( CLOBBER ) then continue

       ;; Order structure
       ordr_str   = mike_getfil('ordr_str', setup, SIDE=qq)

       ;; Process an arc if necessary!
       if ordr_str[0].arc_m[0] EQ 0. then begin
           ;; Find the arc
           print, 'mike_slitflat: Processing an Arc.  Wish me luck!'
           arc = where(mike.type EQ 'ARC' AND mike.flg_anly NE 0 AND $
                       mike.setup EQ setup AND mike.side EQ qq, narc)
           if narc EQ 0 then begin
              print, "mike_slitflat: no arcs to be found"
              continue
           endif

           
           for iarc=0,n_elements(arc) do begin
             fil_nm = mike_getfil('arc_trc', setup, $
                 subfil=mike[arc[iarc]].img_root, /name, chkf=chkf)
             if chkf EQ 1 then begin
               arc_fil = mike[arc[iarc]].rootpth + mike[arc[iarc]].img_root
               break
             endif
           endfor

           ;; Process
           resolve_routine, 'mike_allarc', /no_recompile, /either

           ;; likely we only need to refit arc trace...

           if NOT keyword_set(clobber) then $
             rslt = mike_fittrcarc_work( arc_fil, setup, qq, CLOBBER=1) $
           else  $
             rslt = mike_allarc_sngl(arc_fil, setup, qq, CLOBBER=1)

           if rslt EQ '' then stop
           ;; Reread order structure and check
           ordr_str   = mike_getfil('ordr_str', setup, SIDE=qq)
           if ordr_str[0].arc_m[0] EQ 0. then stop
       endif
           
       ;; TFLAT
       tflat     = mike_getfil('tflat_fil',setup, SIDE=qq)
       tflativar = mike_getfil('tflat_fil',setup, SIDE=qq, INDX=1)

       print, 'mike_slitflat: NOT Applying MFlat..., already done in mktflat'
       ;;
       ;; MFLAT
       ;; mflat = mike_getfil('mflat_fil', setup, SIDE=qq)
       ;; tflat = tflat / (mflat + (mflat LE 0.5))
       ;; tflativar = tflativar * mflat^2 * (mflat GT 0.5) 

       infinite = where(finite(tflat) EQ 0)
       if infinite[0] NE -1 then begin
         tflat[infinite] = 0
         tflativar[infinite] = 0
       endif

       ;  grow the masked pixels in inverse variance with by 3x3
       tflativar = tflativar * (smooth(1.0*(tflativar EQ 0),3) EQ 0)

       if NOT keyword_set(iresidual_fil) then $
         residual_fil = 'Flats/Residual_'+nm+'_'+c_s+'.fits' $
         else residual_fil = iresidual_fil

       if NOT keyword_set(itflat_spec_fil) then $
         tflat_spec_fil = 'Flats/Tflt_SPEC_'+nm+'_'+c_s+'.fits' $
         else tflat_spec_fil = itflat_spec_fil

       ;; Main call
       profile_struct = x_slitflat(tflat, tflativar, ordr_str, $
                                   chk=chk, residual=residual, qa_str=qa_str, $
                                   prof_image=profile_img, $
                                   detilt=(not keyword_set(nodetilt)), $
                                   tflat_spec_fil=tflat_spec_fil)

       ;; Write to fits
       print, 'mike_slitflat: Writing order w/Profile structure ', profile_fil
       mwrfits, profile_struct, profile_fil, /create
       
       ;; Write to fits
       print, 'mike_slitflat: Writing Profile Image structure ', proimg_fil
       mwrfits, profile_img, proimg_fil, /create
       
       if keyword_set(residual) then begin
           print, 'mike_slitflat: Writing Residual Image', residual_fil
           mwrfits, residual, residual_fil, /create
       endif

       qafil = mike_getfil('qa_slitflat', setup, SIDE=qq)
       x_qaslitprof, qa_str, qafil

   endfor

   ;; All done
   print, 'mike_slitflat: All done'
   return

end


