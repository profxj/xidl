;+ 
; NAME:
; sdsscontistrct__define
;   Version 2.0
;
; PURPOSE:
;    Structure for sdss_fndlin continua in CIV, SiIV, CaII searches
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
;   25-May-2011 Created by KLC
;   30-Jul-2011 New slim version, KLC
;-
pro sdsscontistrct__define

  nlin = 300L                   ; max number of centroids allowed
  npix = 4000L                  ; max possible pixels, SDSS ~ 3800
  nconti = 3                    ; number of continua types
  tmp = {sdsscontistrct, $
         ;; QSO information
         qso_name: ' ', $
         MJD: 0L, $             ; spectroscopic
         PLATE: 0L, $
         FIBER: 0L, $
         RA: 0.d, $
         DEC: 0.d, $
         Rmag: 0., $
         sdss_obs: strarr(10), $ ; spectrum file name(s)
         z_qso: 0d, $
         balflg: 0, $                  ; BAL flag (0: is not BAL; >=1: is BAL)
         npix:0L, $                    ; number of pixels in spectrum
         ipix0:0L, $                   ; beginning pixel for search (redward of Lya)
         cflg:-1, $                    ; continuum flag; see sdss_getcflg()
         ;; Observed wavelengths
         ncent:lonarr(nconti), $         ; number of detected lines
         centroid:fltarr(nlin,nconti), $ ; detected line centroids (observed)
         ;; Measurements of line (different from sdsscivstrct)
         wvlim_orig: fltarr(nlin,2), $ ; wavelength bounds; just for cflg item
         ew_orig: fltarr(nlin), $      ; observed
         sigew_orig: fltarr(nlin),$
         ;; Additional spectrum information
         snr_conv:fltarr(npix,nconti), $ ; Gaussian-convolution (x_findgauss)
         conti:dblarr(npix,nconti), $    ; continuum
         sigconti:dblarr(npix,nconti), $ ; continuum error only
         ;; Completeness information (temp)
         deltaz_orig:fltarr(nconti) $ ; 3.5 sigma wavelength range; see snr_conv
        }
  return
end
