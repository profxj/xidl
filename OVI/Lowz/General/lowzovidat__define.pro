;+ 
; NAME:
; lowzovidat__define
;  V1.1
;
; PURPOSE:
;    Structure summarizing the observations for a given field
;
; CALLING SEQUENCE:
;   tmp = {lowzovidat}
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
;   Oct-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro lowzovidat__define


  strctnm = {lowzovidat, qso: '',$
             qso_ra: '',$
             qso_dec: '',$
             qso_zem: 0.0,$
             qso_vmag: 0.0,$
             qso_uv: 0.0,$
             flg_gal: 0,  $
             gal_fil: ' ', $
             R_limit: 0.0,$
             N_gal: lonarr(4), $
             complete: lonarr(4,2), $
             flg_fuse: 0,  $
             fuse_exp: 0.,  $
             fuse_snr: 0.,  $
             flg_stis: 0,  $
             stis_comm: strarr(4),  $
             flg_ghrs: 0,  $
             ghrs_comm: strarr(4)   $
            }

  return
end
