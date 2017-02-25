;+
; NAME:
;   calculate_resolution
;
; PURPOSE:
; calculate_resolution, s, anamorph=anamorph
;
; CALLING SEQUENCE:
;   long_superbias, filenames, outfile, [ sigrej=, maxiter=, /verbose ]
;                 
;
; INPUTS:
;   s -- Spectral structure which contains all the key quantities
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  
; EXAMPLES:
;
; BUGS:
;    
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   11-Mar-2005  Written by J. Hennawi (UCB), D. Schlegel (LBL)
;-
;------------------------------------------------------------------------------
PRO LONG_CALC_RESLN, s, anamorph = anamorph

   IF NOT keyword_set(anamorph) then anamorph = 1.0
   
   nobj =  n_elements(s)
   
   FOR i = 0, nobj-1L DO BEGIN
       spat_bin = s[i].binning[0]
       spec_bin = s[i].binning[1]
       ;; FWHM resolution as predicted by seeing
       spec_fwhm  = djs_median(s[i].FWHMFIT)/anamorph*spat_bin/spec_bin
       ;; FWHM resolution measured from arc lines
       arc_fwhm = s[i].ARC_FWHM_MED
       ;; Resolution cannot be greater than set by the slit width
       s[i].pix_res = (spec_fwhm <  arc_fwhm)
       splog, 'Resolution implied by seeing is ', 100*spec_fwhm/arc_fwhm, $
              ' percent of arc line width'
   ENDFOR
   
   RETURN
END
       




