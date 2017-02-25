;+ 
; NAME:
; hamspec_tracearc   
;     Version 1.1
;
; PURPOSE:
;    To trace the arc lines in each order (individually) and fit a
;    straight line to each one.  The following steps are taken:
;    1.  Scattered light is removed from the image
;    2.  All significant arc lines are identified (5 sigma)
;    3.  trace_crude is used to trace the lines 
;    4.  trace_crude is reapplied to only those lines which are
;    entirely in the order
;    5.  xy2traceset is used to fit a straight line to each arc line
;    6.  Only the good lines are saved for 2D fitting in a structure
;    which is written to disk
;
; CALLING SEQUENCE:
;  hamspec_tracearc, hamspec, setup, [obj_id, chip], /CLOBBER, INIO=
;
; INPUTS:
;   hamspec     -  MIKE structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [chip]   -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  Output structure describing the fits to the arc lines
;
; OPTIONAL KEYWORDS:
;   INIO      - Initial order (for debugging)
;   /CLOBBER  - Overwrite previous fits
;  ARCFIL=   - Name of the arc file to process (Optional to using
;               setup, chip, etc.)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  hamspec_tracearc, hamspec, setup, obj_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  trace_crude
;  xy2traceset
;  x_echtracearc
;
; REVISION HISTORY:
;   28-Apr-2003 Written by SB
;   04-Feb-2013   Modified from HIRES by JXP
;
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_tracearc, hamspec, setup, INIO=inio, CLOBBER=clobber, $
                    ARCFIL=arcfil
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_tracearc, hamspec, setup, INIO=, /CLOBBER'
      print, '     ARCFIL= [v1.1]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( INIO ) then inio = 0L
  if not keyword_set( SAT ) then sat = 50000.


; Loop on chip
  if not keyword_set(ARCFIL) then begin
     ;; Grab all obj indices
     arcfil = hamspec_getfil('arc_fil', setup, /name)
  endif 

  ;; ORD_STR
  ordr_str = hamspec_getfil('ordr_str', setup)

  
  ;;  Check for outfil
  out_fil = hamspec_getfil('arc_trc', setup, $
                           /name, CHKFIL=chkf)
  ;;  Check for outfil
  if chkf NE 0 and not keyword_set( CLOBBER ) then begin
     print, 'hamspec_tracearc: Arc fit file exists. ' + $
            'Continuing..'
     return
  endif
  ;; QA
  qafil = hamspec_getfil('qa_tracearc', setup, $
                         /name, CHKFIL=chkf)
  rslt = x_echtrcarc(arcfil, ordr_str, out_fil, $
                     CLOBBER=clobber, INIO=inio, QAFIL=qafil)
  if size(rslt,/tname) NE 'STRING' then stop

  print, 'hamspec_tracearc:  All done!'
  return
end

