;+ 
; NAME:
; emissstrct__define
;   Version 1.1
;
; PURPOSE:
;  Structure for a simple emission line. 
;
; CALLING SEQUENCE:
;   tmp = {emissstrct}
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
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;  Written by JXP
;-
;------------------------------------------------------------------------------
pro emissstrct__define

;  This routine defines the line list structure

  tmp = {emissstrct, $
         ion: ' ', $
         wrest: 0.d, $
         f: 0.d, $
         gamma: 0., $
         wcen: 0.d, $    ;; Observed wavelength centroid
         sig_wc: 0.d, $
         vsigma: 0.d, $
         sig_vsig: 0.d, $
         A_gauss: 0., $  ;; Height of Gaussian
         sig_A: 0., $
         set: 0, $              ; Groups abs lines together
         flg_flux: 0, $         ; Flag  (0=unmeasured, 1=Gaussian, 2=Boxcar)  [Negative for upper limit]
         flux: 0.d, $               ; Colm (log) 
         fsig: 0.d, $
         flg_c: 0, $           ; Flag for continuum mode (0=constant; 1=linear)
         c_reg: dblarr(2,2), $  ; Continuum regions
         f_reg: dblarr(2), $  ; Fit region
         EW: 0., $
         EWsig: 0., $
         zem: 0.d,  $          ; Redshift
         zsig: 0.d  $           ; Error in redshift
        }

end
  
         
