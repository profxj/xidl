;+ 
; NAME:
; lowzovidat__define
;  V1.1
;
; PURPOSE:
;    Given a list of DLA base files, fill up the structure ;
; CALLING SEQUENCE:
;   
;   lowzovi_prsdat, stucture, filename
;
; INPUTS:
;
; RETURNS:
;   structure      - IDL structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  LIST - File
;  ION - Input ionic column densities
;  NOELM - Supress inputting Elemental values
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lowzovi_prsdat, struct, '/u/xavier/DLA/Lists/tot_dla.lst'
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   26-Aug-2003 Written by JXP
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
