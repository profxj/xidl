;+ 
; NAME:
; sdssmcstrct__define
;   Version 1.0
;
; PURPOSE: 
;    Structure for completeness test intermediate information
;    per Monte Carlo realization per sightline
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
;   21-Sep-2011 Created by KLC
;-

pro sdssmcstrct__define
  ncent = 300L                  ; match sdsscontistrct
  nlin = 10                     ; max number of input lines 
;  npix = 4000L                  ; max possible pixels, SDSS ~ 3800
  ncomp = 50                    ; max allowed number of components
  tmp = {sdssmcstrct, $   
         ;; QSO information
         qso_name: ' ', $
         MJD: 0L, $             ; spectroscopic
         PLATE: 0L, $
         FIBER: 0L, $
         RA: 0.d, $
         DEC: 0.d, $
         Rmag: 0., $
         z_qso: 0d, $
         balflg: 0, $           ; BAL flag (0: is not BAL; >=1: is BAL)
         cflg:-1, $             ; continuum flag; see sdss_getcflg
         ;; Scrubbed centroids
         cent_mask:intarr(ncent), $
         snr:0., $              ; of new cleaned spectrum
;         pix_mask:intarr(npix), $ ; 0: good; 1: bad
         ;; Simulated doublets
         ;; see use in sdss_genprof_cpvp2mcstrct()
         nsys:0, $              ; number of systems in arrays
         ncomp:intarr(nlin), $
         id_sys:lonarr(nlin), $ ; unique ID number
         wrest:fltarr(nlin,2), $
         wvlim_input:fltarr(nlin,2,2), $ 
         zabs_input:fltarr(nlin,ncomp+1), $ 
         ncolm_input:fltarr(nlin,ncomp+1), $ ; Input number
         b_input:fltarr(nlin,ncomp+1), $ ; Doppler param (km/s)
         ew_input:fltarr(nlin,2), $ 
         ;; Recovered params
         flg_rec:intarr(nlin), $ ; 0 = not rec; 1 = recovered
         lsnr:fltarr(nlin,2), $ ; LSNR of lines (sdss_fndlin, sdss_fndciv)
         wvlim_rec:fltarr(nlin,2,2), $ 
         zabs_rec:fltarr(nlin,2), $ 
         sigzabs_rec:fltarr(nlin,2), $ 
         ncolm_rec:fltarr(nlin,2), $ ; AODM
         signcolm_rec:fltarr(nlin,2), $ 
         ncolmflg_rec:intarr(nlin,2), $
         ew_rec:fltarr(nlin,2), $ 
         sigew_rec:fltarr(nlin,2) $ 
         }

  return

end 
