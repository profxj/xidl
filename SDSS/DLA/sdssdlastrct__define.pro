;+ 
; NAME:
; sdssdlastrct__define
;
; PURPOSE:
;    Structure for SDSS DLA searches
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
;   27-Feb-2004 Written by SHF
;   14-Oct-2004 Modified by JXP
;-
;------------------------------------------------------------------------------
pro sdssdlastrct__define

tmp={sdssdlastrct, $
     qso_name: ' ', $
     plate: 0L, $
     fiber: 0L, $
     mjd: 0L, $
     RA: 0.d, $
     DEC: 0.d, $
     Rmag: 0., $
     sdss_obs: strarr(10), $
     NHI: 0., $
     z_qso: 0., $
     zabs: 0.d, $
     conti: 0., $
     flg_mtl: 0, $
     flg_NHI: 0, $
     quality: 0. $
     }
return
end
     
