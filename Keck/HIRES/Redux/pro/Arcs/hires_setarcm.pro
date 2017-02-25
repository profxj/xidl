;+ 
; NAME:
; hires_setarcm
;     Version 1.1
;
; PURPOSE:
;  Run through the necessary steps to set the arc_m values in the
;  order structure.  This seperate routine was created to avoid the
;  high level of confusion occuring in hires_allarc.
;
; CALLING SEQUENCE:
;   
;   hires_setarcm, arc_fil, setup, chip
;
; INPUTS:
;   arc_fil -  Raw arc file used to set the arc_m values
;   setup   -  Integer defining setup
;   chip    -  Blue (1) or Red (2)
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
;   hires_setarcm, 'Raw/hires_mb0032', 1, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_procarc
;  hires_alignarc
;  hires_fitarc
;  hires_fit2darc
;  hires_tracearc
;  hires_fittrcarc
;
; REVISION HISTORY:
;   18-Aug-2004 Written by JXP
;
;  Usage:
;-
;

pro hires_setarcm, hires, fil, setup, chip, CHK=chk, XYOFF=xyoff, $
                  PINTER=pinter, BAD_ANLY=bad_anly, ARC_FIL=arc_fil, $
                  FITPRM=fitprm, EXTEN=exten

;

  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'hires_setarcm, hires, fil, setup, chip [v1.1]'
      return
  endif 

  ;; Run through the steps without aligning
  ;; Process
  rslt = hires_procarc_sngl( fil, setup, chip, exten, /CLOBBER)
  arc_fil = rslt
  ;; Fit 1D
  rslt = hires_fitarc_work( hires, setup, ARCFIL=arc_fil, $
                            /CLOBBER, PINTER=pinter)
  ;; Fit 2D
  rslt = hires_fit2darc_work( hires, setup, ARCFIL=arc_fil, /CLOBBER)

  ;; Trace
  rslt = hires_tracearc_work( arc_fil, setup, ARCFIL=arc_fil, /CLOBBER)

  ;; Fit Trace
  rslt = hires_fittrcarc_work( arc_fil, setup, ARCFIL=arc_fil, $
                               /CLOBBER, /ORDRCLOB)

  ;; Align
;  ordr_str = hires_getfil('ordr_str', setup, SIDE=chip)
;  arc = xmrdfits(arc_fil, 0, /silent)
;  arci = xmrdfits(arc_fil, 1, /silent)

;  if keyword_set(chk) then window, 0, title='hires_arcalign'
;  xyoff = hires_arcalign_work(arc, arci, ordr_str, CHK=keyword_set(chk), FITPRM=fitprm)

  ;; Rerun steps with alignment in place
  ;; Fit 1D
;  rslt = hires_fitarc_work( arc_fil, setup, chip, /CLOBBER, PINTER=pinter, $
;                         SHFTPRM=fitprm)
  ;; Fit 2D
;  rslt = hires_fit2darc_work( arc_fil, setup, chip, /CLOBBER)
  ;; Trace
;  rslt = hires_tracearc_work( arc_fil, setup, chip, /CLOBBER, SHFTPRM=fitprm)
  ;; Fit Trace
;  rslt = hires_fittrcarc_work( arc_fil, setup, chip, /CLOBBER, /ORDRCLOB)

  ;; Re-Align
;  ordr_str = hires_getfil('ordr_str', setup, SIDE=chip)
;  if keyword_set(chk) then window, 0, title='hires_arcalign'
;  xyoff = hires_arcalign_work(arc, arci, ordr_str, CHK=keyword_set(chk), FITPRM=fitprm)

  print, 'hires_setarcm: All done!'

  return
end

