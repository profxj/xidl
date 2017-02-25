;+ 
; NAME:
; mthr_thruput
;    Version 1.1
;
; PURPOSE:
;    Estimates the throughput of MTHR using theoretical estimates
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  OUTDIR=  -- Name of output directory
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  showfits
;  querydss
;
; REVISION HISTORY:
;   27-Oct-2005 Written by JXP based on HIRES S2N code
;-
;------------------------------------------------------------------------------
function mthr_thruput, wave

;
;..............................................................................
;
  if not keyword_set(icd) then icd = 0 
 
  xwave = [3200, 3500, 4000, 4500, 5000, 6000, 6500., 7000, 8000, 9000., 9500]

  ;;  MTHR estimates, November 2005
  xthru = [8.87604E-08,	0.101181256, 0.146305043,0.14812337,0.18782876,	$
           0.217278499,0.213093159,0.205556217,0.235631684,0.169817629, $
           0.104468766]

  thru = interpol(xthru, xwave, wave)
 
  return, thru
end

















































































































































