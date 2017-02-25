;+ 
; NAME:
; uves_slitflat
;     Version 1.1
;
; PURPOSE:
;    Stores slit profile and gradient along each order from twilight flats.
;    This routine is critical for performing ideal sky subtraction,
;    especially given the short slit length.  The following steps are
;    performed in uves_slitflat_work:
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
;  uves_slitflat, uves, setup, side, [/chk, /clobber]
;
; INPUTS:
;   uves     -  MIKE structure
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
;   uves_slitflat, uves, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_ordermask
;  uves_fitgap
;  uves_qw
;  uves_slitflat_work
;
; REVISION HISTORY:
;   ??-??-2003 Written by SB
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_slitflat, uves, setup, side, chk=chk, $
            tflat_fil=tflat_fil, profile_fil=iprofile_fil, $
            residual_fil=iresidual_fil, CLOBBER=clobber, NODETILT=nodetilt, $
            tflat_spec_fil=itflt_spec_fil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'uves_slitflat, uves, setup, [side], /clobber, /chk, /NODETILT' 
      print, '    RESIDUAL_FIL=, TFLAT_FIL=, PROFILE_FIL=  [v1.1]'
      return
  endif

   if not keyword_set( SIDE ) then side = [1L]
   if not keyword_set( chk  ) then chk=0

; Setup
   if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 

;  Loop on side
   for ii=0L,n_elements(side)-1 do begin
       qq = side[ii]
       ;; SIDE
       if qq EQ 1 then begin
           print, 'uves_slitflat: Profiling BLUE trace flat'
           nm = 'B'
       endif else begin
           print, 'uves_slitflat: Profiling RED trace flat'
           nm = 'R' 
       endelse

      ;; Set wcen
      flt = where(uves.side EQ qq AND uves.flg_anly NE 0 AND $
                  strtrim(uves.type,2) EQ 'TFLT' and $
                  uves.setup EQ setup, nflt)
      wcen = uves[flt[0]].xdangl

      ;;
      if NOT keyword_set(iprofile_fil) then $
        profile_fil = uves_getfil('ordr_str', setup, WCEN=wcen, /name) $
      else profile_fil = iprofile_fil
      
      if NOT keyword_set(iproimg_fil) then $
        proimg_fil = 'Flats/Profile_'+nm+'_'+c_s+'.fits' $ 
      else proimg_fil = iproimg_fil
      if x_chkfil(proimg_fil+'*',/silent) NE 0 $
        and not keyword_set( CLOBBER ) then continue

      ;; Order structure
      ordr_str = uves_getfil('ordr_str', setup, WCEN=wcen)
      
      ;; Process an arc if necessary!
      if ordr_str[0].arc_m[0] EQ 0. then begin
          stop
          ;; Find the arc
          print, 'uves_slitflat: Processing an Arc.  Wish me luck!'
          arc = where(uves.type EQ 'ARC' AND uves.flg_anly NE 0 AND $
                      uves.setup EQ setup AND uves.side EQ qq, narc)
          if narc EQ 0 then begin
              print, "uves_slitflat: no arcs to be found"
              continue
          endif
          
           
           for iarc=0,n_elements(arc) do begin
             fil_nm = uves_getfil('arc_trc', setup, $
                 subfil=uves[arc[iarc]].img_root, /name, chkf=chkf)
             if chkf EQ 1 then begin
               arc_fil = uves[arc[iarc]].rootpth + uves[arc[iarc]].img_root
               break
             endif
           endfor

           ;; Process
           resolve_routine, 'uves_allarc', /no_recompile, /either

           ;; likely we only need to refit arc trace...

;           if NOT keyword_set(clobber) then $
;             rslt = uves_fittrcarc_work( arc_fil, setup, qq, CLOBBER=1) $
;           else  $
;             rslt = uves_allarc_sngl(arc_fil, setup, qq, CLOBBER=1)

           if rslt EQ '' then stop
           ;; Reread order structure and check
           ordr_str   = uves_getfil('ordr_str', setup, SIDE=qq)
           if ordr_str[0].arc_m[0] EQ 0. then stop
       endif
           
       ;; TFLAT
       tflat = uves_getfil('qtz_fil', setup, WCEN=wcen)
       tflativar = uves_getfil('qtz_fil', setup, WCEN=wcen, indx=1)

       print, 'uves_slitflat: NOT Applying MFlat..., already done in mktflat'
       ;;
       ;; MFLAT
       ;; mflat = uves_getfil('mflat_fil', setup, SIDE=qq)
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
       print, 'uves_slitflat: Writing order w/Profile structure ', profile_fil
       mwrfits, profile_struct, profile_fil, /create
       
       ;; Write to fits
       print, 'uves_slitflat: Writing Profile Image structure ', proimg_fil
       mwrfits, profile_img, proimg_fil, /create
       
       if keyword_set(residual) then begin
           print, 'uves_slitflat: Writing Residual Image', residual_fil
           mwrfits, residual, residual_fil, /create
       endif

       qafil = uves_getfil('qa_slitflat', setup, WCEN=wcen)
       x_qaslitprof, qa_str, qafil

   endfor

   ;; All done
   print, 'uves_slitflat: All done'
   return

end
