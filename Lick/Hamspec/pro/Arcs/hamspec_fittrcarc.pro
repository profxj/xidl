;+ 
; NAME:
; hamspec_fittrcarc
;     Version 1.1
;
; PURPOSE:
;   To fit the slope of the arc lines as a function of order number
;   and y position on the CCD.  This information is then used to
;   construct a 2D wavelength image.  The fitting routine is the usual
;   least-squares algorithm with two rounds of rejection.
;
; CALLING SEQUENCE:
; hamspec_fittrcarc, hamspec, setup, obj_id, chip, /CHK, /CLOBBER, $
;                   ARCFIL=arcfil
;
; INPUTS:
;   hamspec    -  HIRES structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [chip]   -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  Fits file with the coefficients of the 2D fit.  Filename like
;  'Arcs/TRC/Arc_B0539_F.fits' 
;
; OPTIONAL KEYWORDS:
;  /CHK  -- Plots residuals
;  /CLOBBER -- Overwrite previous solution
;  ARCFIL=   - Name of the arc file to process (Optional to using
;               setup, chip, etc.)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Feb-2005 Written by JXP
;   04-Feb-2013   Modified from HIRES by JXP
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_fittrcarc, hamspec, setup, CHK=chk, CLOBBER=clobber, $
                     ARCFIL=arcfil
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_fittrcarc, hamspec, setup, /CHK, /CLOBB, ARCFIL= [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set(ARCFIL) then begin
     arcfil = hamspec_getfil('arc_fil', setup, /name)
  endif

  ;; ORD_STR
  ordr_str = hamspec_getfil('ordr_str', setup, fil_nm=ordr_fil)

  ;;  Check for outfil
  out_fil = hamspec_getfil('arc_fittrc', setup, /name, CHKFIL=chkf)  
  ;;  Check for outfil
  if chkf NE 0 and not keyword_set( CLOBBER ) then begin
     print, 'hamspec_fittrcarc: Arc fit file exists. ' + $
            'Continuing..'
     return
  endif

  ;; TRC_FIL
  trc_fil = hamspec_getfil('arc_trc', setup, $
                           /name, CHKFIL=chkf) 
  ;; QA
  qafil = hamspec_getfil('qa_fittrcarc', setup, $
                         /name, CHKFIL=chkf)
  
  ;; Main driver
  rslt = x_fittrcarc(arcfil, trc_fil, ordr_str, out_fil, qafil, $
                     CHK=chk, CLOBBER=clobber, $
                     ORDR_FIL=ordr_fil) 
  if size(rslt,/tname) NE 'STRING' then stop

  print, 'hamspec_fittrcarc: All done!'

  return
end
