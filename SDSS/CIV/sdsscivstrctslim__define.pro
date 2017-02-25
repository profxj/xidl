;+ 
; NAME:
; sdsscivstrctslim__define
;   Version 2.0
;
; PURPOSE:
;    Structure for CIV search
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
;   24-May-2011 Adopted from sdssmgiistrct__define, KLC
;   13-Nov-2011 Adopted from sdsscivstrct__define, KLC
;-
pro sdsscivstrctslim__define
  nlin = 2
  tmp={sdsscivstrctslim, $
       ;; QSO information (RA, DEC, Rmag in sdsscontistrct)
       qso_name: ' ', $
       MJD: 0L, $               ; spectroscopic
       PLATE: 0L, $
       FIBER: 0L, $
       sdss_obs: '', $  ; spectrum file name(s)
       z_qso: 0d, $
       balflg: 0, $             ; BAL flag (0: is not BAL; 1: is BAL)
       ;; Absorber information (different from sdsscontistrct)
       wrest: dblarr(nlin),  $        ; rest wavelength
       wvlim_final: dblarr(nlin,2), $ ; wavelength bounds [0,*] is 1548, [1,*] is 1550
       wvlim_orig: dblarr(nlin,2), $  ; from sdss_ewciv
       ;; Redshift
       zabs_final: dblarr(nlin), $ ; used for publication
       sigzabs_final: dblarr(nlin), $
       zabs_orig: dblarr(nlin), $ ; original from sdss_civsearch
       gwidth: fltarr(nlin), $    ; hold LSNR cut per line
       ;; REST Equivalent widths (_orig and _tau in sdsscontistrct are
       ;; observed and just different from here)
       ew_final: fltarr(nlin), $ ; [1548, 1550]
       sigew_final: fltarr(nlin),$
       ew_orig: fltarr(nlin), $ ; original from sdss_civsearch [boxcar]
       ewflg: intarr(nlin), $   ; Limit; see sdss_getlimflg()
       ;; log column density 
       ncolm_final: fltarr(nlin), $
       signcolm_final: fltarr(nlin),$
       ncolm_orig: fltarr(nlin), $ ; AODM from sdss_civsearch [boxcar]
       ncolmflg: intarr(nlin), $ ; Limit; see sdss_getlimflg()
       ;; Misc.
       rating: replicate(-1,10), $ ; from sdss_chkciv for system; rating[9] = blend flag (1 = yes)
       abslin_fil:'', $            ; ABSLIN/ file from sdss_fndlin; see sdsscontistrct
       cflg: -1, $                 ; continuum flag; see sdss_getcflg()
       cmplt_fil:'', $             ; completeness file name
       notes:'' $                  ; comments
      }
  return
end
     
