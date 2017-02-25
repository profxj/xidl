;+ 
; NAME:
; hamspec_mkaimg
;     Version 1.1
;
; PURPOSE:
;   Given the 2D solution for the slope of the lines as a function of
;   position this code creates a wavelength image (i.e. assigns a
;   unique wavelength to each pixel in each order).  The xoffset is 
;   input in order to properly determine the edges of each order.  A
;   simple spline interpolation is used to determine the values.
;
;
; CALLING SEQUENCE:
;  hamspec_mkaimg, hamspec, setup, [obj_id, chip], /CHK, /CLOBBER, ARCFIL=
;
; INPUTS:
;   hamspec     -  MIKE structure
;   setup    -  Integer defining setup
;   [obj_id]   -  Object identifier
;   [chip]   -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  2D wavelength image with name like 'Arcs/Arc_mb0439I.fits'
;
; OPTIONAL KEYWORDS:
;  /CLOBBER  - Overwrite previous image
;  /CHK      - Display the final image
;  ARCFIL=   - Name of the arc file to process (Optional to using
;               setup, chip, etc.)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_mkaimg, hamspec, setup, obj_id, chip
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-May-2003 Written by SB
;   04-Feb-2013   Modified from HIRES by JXP
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hamspec_mkaimg, hamspec, setup, CHK=chk, CLOBBER=clobber, $
                  ARCFIL=arcfil, BAD_ANLY=bad_anly, _EXTRA=extra

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_mkaimg, hamspec, setup, /CHK, /CLOBBER, ARCFIL= [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set(ARCFIL) then begin
     arcfil = hamspec_getfil('arc_fil', setup, /name)
  endif

  ;; ORD_STR
  ordr_str = hamspec_getfil('ordr_str', setup, fil_nm=ordr_fil)

  ;;  Check for outfil
  out_fil = hamspec_getfil('arc_mkaimg', setup, $
                           /name, CHKFIL=chkf) 
  ;;  Check for outfil
  if chkf NE 0 and not keyword_set( CLOBBER ) then begin
     print, 'hamspec_fittrcarc: Arc fit file exists. ' + $
            'Continuing..'
     return
  endif

  arc2d_fil = hamspec_getfil('arc_2Dfit', setup, $
                             /name, CHKFIL=chkf )
  fil_fittrc = hamspec_getfil('arc_fittrc', setup, $
                                 /name, CHKFIL=chkf)
;          stop ;; Odds are that arc_xyoff is not set!
  rslt = x_mkaimg(arcfil, ordr_str, arc2d_fil, fil_fittrc, $
                  out_fil, CHK=chk, CLOBBER=clobber, $
;                              SHFTPRM=hamspec[idx].arc_xyoff, $
                  BAD_ANLY=bad_anly, _EXTRA=extra) 

  if bad_anly GT 0 then stop
;             remove = where(strpos(hamspec.img_root, arc_roots[idx]) GT -1)
;             if remove[0] NE -1 then hamspec[remove].flg_anly = 0
;          endif

  if size(rslt,/tname) NE 'STRING' then stop
  print, 'hamspec_mkaimg: All done!'

  return
end
