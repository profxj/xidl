;+ 
; NAME:
; mike_setarcm
;     Version 1.1
;
; PURPOSE:
;  Run through the necessary steps to set the arc_m values in the
;  order structure.  This seperate routine was created to avoid the
;  high level of confusion occuring in mike_allarc.
;
; CALLING SEQUENCE:
;   
;   mike_setarcm, arc_fil, setup, side
;
; INPUTS:
;   arc_fil -  Raw arc file used to set the arc_m values
;   setup   -  Integer defining setup
;   side    -  Blue (1) or Red (2)
;
; RETURNS:
;
; OUTPUTS:
;  A series of files related to Arc calibration
;
; OPTIONAL KEYWORDS:
;  /PINTER   - Perform fit for pre-identified lines
;
; OPTIONAL OUTPUTS:
;  ARC_FIL=  - Name of processed arc
;  XYOFF=    - Offset fit info
;
; ADDITIONAL KEYWORDS TO MIKE_ALLARC_SNGL:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_setarcm, 'Raw/mike_mb0032', 1, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_procarc
;  mike_alignarc
;  mike_fitarc
;  mike_fit2darc
;  mike_tracearc
;  mike_fittrcarc
;
; REVISION HISTORY:
;   18-Aug-2004 Written by JXP
;
;  Usage:
;-
;

pro mike_setarcm, fil, setup, side, CHK=chk, XYOFF=xyoff, $
                  PINTER=pinter, BAD_ANLY=bad_anly, ARC_FIL=arc_fil, $
                  FITPRM=fitprm, _EXTRA=extra

;

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'mike_setarcm, fil, setup, side [v1.1]'
      return
  endif 

  ;; Run through the steps without aligning
  ;; Process
  rslt = mike_procarc_sngl( fil, setup, side, _EXTRA=extra, /CLOBBER)
  arc_fil = rslt
  ;; Fit 1D
  rslt = mike_fitarc_work( arc_fil, setup, side, /CLOBBER, PINTER=pinter)
  ;; Fit 2D
  rslt = mike_fit2darc_work( arc_fil, setup, side, /CLOBBER)
  ;; Trace
  rslt = mike_tracearc_work( arc_fil, setup, side, /CLOBBER)
  ;; Fit Trace
  rslt = mike_fittrcarc_work( arc_fil, setup, side, /CLOBBER, /ORDRCLOB)

  ;; Align
  ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
  arc = xmrdfits(arc_fil, 0, /silent)
  arci = xmrdfits(arc_fil, 1, /silent)

  if keyword_set(chk) then window, 0, title='mike_arcalign'
  xyoff = mike_arcalign_work(arc, arci, ordr_str, CHK=keyword_set(chk), FITPRM=fitprm)

  ;; Rerun steps with alignment in place
  ;; Fit 1D
  rslt = mike_fitarc_work( arc_fil, setup, side, /CLOBBER, PINTER=pinter, $
                         SHFTPRM=fitprm)
  ;; Fit 2D
  rslt = mike_fit2darc_work( arc_fil, setup, side, /CLOBBER)
  ;; Trace
  rslt = mike_tracearc_work( arc_fil, setup, side, /CLOBBER, SHFTPRM=fitprm)
  ;; Fit Trace
  rslt = mike_fittrcarc_work( arc_fil, setup, side, /CLOBBER, /ORDRCLOB)

  ;; Re-Align
  ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
  if keyword_set(chk) then window, 0, title='mike_arcalign'
  xyoff = mike_arcalign_work(arc, arci, ordr_str, CHK=keyword_set(chk), FITPRM=fitprm)

  print, 'mike_setarcm: All done!'

  return
end

