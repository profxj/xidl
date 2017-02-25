;+ 
; NAME:
; sdssmgiistrct__define
;
; PURPOSE:
;    Structure for MgII search
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
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   27-Feb-2004 Written by GEP
;-
pro sdssmgiistrct__define

tmp={sdssmgiistrct, $
     qso_name: ' ', $
     PLATE: 0L, $
     FIBER: 0L, $
     RA: 0.d, $
     DEC: 0.d, $
     Rmag: 0., $
     sdss_obs: strarr(10), $
     z_qso: 0., $
     zabs: 0.d, $
     gwidth: 0., $
     ew: fltarr(50), $
     sigew: fltarr(50),$
     wrest: dblarr(50),  $
     ewflg: 0 $
     }
return
end
     
