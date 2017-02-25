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
pro sdssllsstrct__define

tmp={sdssllsstrct, $
     plate: 0L, $
     mjd: 0L, $
     datfil: '', $
     gmag: 0., $
     umag: 0., $
     uerr: 0., $
     fiber: 0L, $
     zem: 0.0, $
     zlls: 0.0, $
     taulls: 0.0, $
     sigtau: 0.0, $
     blls: 0.0, $
     llsflg: 0, $
     flg_NHI: -1, $
     flg_LLS: -1, $
     flg_tweak: intarr(3), $
     flg_extra: 0, $
     conti_scale: 0., $
     zstart: 0.0, $
     zend: 0.0, $
     svchi: fltarr(100), $
     svz: fltarr(100), $
     uminusg: 0.0} 

return
end
     
