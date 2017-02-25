;+ 
; NAME:
; wfccd_allarcsol
;    Version 1.0
;
; PURPOSE:
;    Solves arc solutions for a given mask
;      Designed to do only 1 at a time
;
; CALLING SEQUENCE:
;   
;   wfccd_allarcsol, wfccd
;
; INPUTS:
;   wfccd     - WFCCD structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL INPUTS:
;   w0off      - half-range of lambda zeropoint to check on initial guess
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_allarcsol, wfstrct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


pro wfccd_allarcsol, wfccd, mask_id, NOFITS=nofits, CLOBBER=clobber, $
                     W0off=w0off, W2OFF=w2off

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_allarcsol, wfccd, mask_id, /NOFITS [v1.0]'
    return
  endif 

;  Optional Keywords

;  Find OBJ frames
  obj = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nobj)
  arcfil = wfccd[obj].arc_fil

;  LOOP ON Unique Arc Files
  for q=0L,nobj-1 do $
    wfccd_arcsol, wfccd, obj[q], NOFITS=nofits, CLOBBER=clobber, W0OFF=w0off, $
    W2OFF=w2off

  return
end
