;+ 
; NAME:
; hamspec_slitflat
;     Version 1.1
;
; PURPOSE:
;    Stores slit profile and gradient along each order from twilight flats.
;    This routine is critical for performing ideal sky subtraction,
;    especially given the short slit length.  The following steps are
;    performed in hamspec_slitflat_work:
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
;  hamspec_slitflat, hamspec, setup, chip, [/chk, /clobber]
;
; INPUTS:
;   hamspec    -  HIRES structure
;   setup    -  Setup identifier 
;   [chip]   -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
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
;  /NODETILT   - Do NOT remove a linear tilt from the Cross-section fit
;
;  NXBKT  -  Number of x breakpoints for scattered light fit 
;  NYBKT  -  Number of x breakpoints for scattered light fit 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_slitflat, hamspec, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;  hamspec_ordermask
;  hamspec_fitgap
;  hamspec_qw
;  hamspec_slitflat_work
;
; REVISION HISTORY:
;   ??-??-2003 Written by SB
;   04-Feb-2013   Modified from HIRES by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_slitflat, hamspec, allsetup, chk=chk, $
            tflat_fil=tflat_fil, profile_fil=iprofile_fil, $
            residual_fil=iresidual_fil, CLOBBER=clobber, NODETILT=nodetilt, $
            tflat_spec_fil=itflt_spec_fil, ZERO_SLOPE=zero_slope

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_slitflat, hamspec, setup, [chip], /clobber, /chk, /NODETILT' 
      print, '    RESIDUAL_FIL=, TFLAT_FIL=, PROFILE_FIL=, /ZERO_SLOPE  [v1.1]'
      return
  endif

   if not keyword_set( chk  ) then chk=0

   for kk=0L,n_elements(allsetup)-1 do begin
       setup = allsetup[kk]
       ;; Setup
       if setup GE 10 then c_s = strtrim(setup,2) $
       else c_s = '0'+strtrim(setup,2) 

       if NOT keyword_set(iprofile_fil) then $
          profile_fil = hamspec_getfil('ordr_str', setup, /name) $
       else profile_fil = iprofile_fil
           
       if NOT keyword_set(iproimg_fil) then $
          proimg_fil = 'Flats/Profile_'+c_s+'.fits' $ 
       else proimg_fil = iproimg_fil
       if x_chkfil(proimg_fil+'*',/silent) NE 0 $
          and not keyword_set( CLOBBER ) then continue
       
       ;; Order structure
       ordr_str   = hamspec_getfil('ordr_str', setup, fil_nm=ordr_fil)
       
       ;; Process an arc if necessary!
       if ordr_str[0].arc_m[0] EQ 0. and not keyword_set(ZERO_SLOPE) then begin
          ;; Do not do this anymore
          print, 'hamspec_slitflat: Process an arc first and then return'
          print, '   to this step.  If you have already processed an arc,'
          print, '   then reprocess it because you have managed to overwrite'
          print, '   the order structure. [Use /CLOBBER]'
          print, '   Returning....'
          return
          
          ;; Find the arc
          print, 'hamspec_slitflat: Processing an Arc.  Wish me luck!'
          arc = where(hamspec.type EQ 'ARC' AND hamspec.flg_anly NE 0 AND $
                      hamspec.setup EQ setup AND hamspec.chip EQ qq, narc)
          if narc EQ 0 then begin
             print, "hamspec_slitflat: no arcs to be found"
             continue
          endif
          
          
          for iarc=0,n_elements(arc) do begin
             pos = strpos(hamspec[arc[iarc]].arc_fil, '.fits')
             frame = long(strmid(hamspec[arc[iarc]].arc_fil,pos-4,4))
             fil_nm = hamspec_getfil('arc_trc', setup, CHIP=qq, $
                                     /name, CHKFIL=chkf,  FRAME=frame)
             if chkf EQ 1 then begin
                arc_fil = hamspec[arc[iarc]].rootpth + hamspec[arc[iarc]].img_root
                break
             endif
          endfor
               
          ;; likely we only need to refit arc trace...
          if NOT keyword_set(clobber) then begin
             trc_fil = hamspec_getfil('arc_trc', setup, $
                                      /name, CHKFIL=chkf)
             qafil = hamspec_getfil('qa_fittrcarc', setup, $
                                    /name, CHKFIL=chkf)
             out_fil = hamspec_getfil('arc_fittrc', setup, $
                                      /name, CHKFIL=chkf)
             rslt = x_fittrcarc(arc_fil, trc_fil, ordr_str, out_fil, qafil, $
                                CHK=chk, CLOBBER=clobber, $
                                ORDR_FIL=ordr_fil) 
          endif else  stop
          if rslt EQ '' then stop
          ;; Reread order structure and check
          ordr_str   = hamspec_getfil('ordr_str', setup)
          if ordr_str[0].arc_m[0] EQ 0. then stop
       endif
           
       ;; TFLAT
       tflat     = hamspec_getfil('qtz_fil',setup)
       tflativar = hamspec_getfil('qtz_fil',setup, INDX=1)
       scatt_img = hamspec_getfil('nqtz_fil',setup, INDX=2, /name,$
                                   CHKFIL=chkf)
       if CHKF EQ 0 then begin
          ;hamspec_rdxlog, 'flat_log', $
          ;               'hamspec_slitflat: Subtracting scattered light'
          sz = size(tflat,/dimensions)
          ;; Calculated scattered light image
          msk = x_ordermask(sz[0], sz[1], ordr_str, trim=1)
          nxbkpt=5L 
          nybkpt=5L
          scatt_img = x_fitgap(tflat, tflativar, msk, $
                               nxbkpt=nxbkpt, nybkpt=nybkpt)
       endif
       
       ;; Subtract scattered light image
       tflat = tflat - scatt_img
       
       ;; MFLAT
       
       infinite = where(finite(tflat) EQ 0)
       if infinite[0] NE -1 then begin
          tflat[infinite] = 0
          tflativar[infinite] = 0
       endif
       
       ;;  grow the masked pixels in inverse variance with by 3x3
       tflativar = tflativar * (smooth(1.0*(tflativar EQ 0),3) EQ 0)
       
       if NOT keyword_set(iresidual_fil) then $
          residual_fil = 'Flats/Residual_'+c_s+'.fits' $
       else residual_fil = iresidual_fil
       
       if NOT keyword_set(itflat_spec_fil) then $
          tflat_spec_fil = 'Flats/Tflt_SPEC_'+c_s+'.fits' $
       else tflat_spec_fil = itflat_spec_fil
       
       ;; Set breakpoint spacing  (bigger means more breakpoints)
       bluscl = 3.
       
       ;; Main call
       profile_struct = x_slitflat(tflat, tflativar, ordr_str, $
                                   chk=chk, residual=residual, $
                                   qa_str=qa_str, $
                                   prof_image=profile_img, $
                                   detilt=(not keyword_set(nodetilt)), $
                                   tflat_spec_fil=tflat_spec_fil,$
                                   /NOSCATT)
;                                       bluscl = bluscl, /NOSCATT, $
;                                       SCLLIMIT=1000.)
       
       ;; Write to fits
       print, 'hamspec_slitflat: Writing order w/Profile structure ', profile_fil
       mwrfits, profile_struct, profile_fil, /create
       
       ;; Write to fits
       print, 'hamspec_slitflat: Writing Profile Image structure ', proimg_fil
       mwrfits, profile_img, proimg_fil, /create
       
       if keyword_set(residual) then begin
          print, 'hamspec_slitflat: Writing Residual Image', residual_fil
          mwrfits, residual, residual_fil, /create
       endif
       
       qafil = hamspec_getfil('qa_slitflat', setup)
       x_qaslitprof, qa_str, qafil
    endfor
           
   ;; All done
   print, 'hamspec_slitflat: All done'
   return

end


