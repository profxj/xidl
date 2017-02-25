;+ 
; NAME:
; sdss_llsinit
;    Version 1.0
;
; PURPOSE:
;    Produces a LLS structure of SDSS LLS 
;
; CALLING SEQUENCE:
;  sdss_llsstrct, lls
;
; INPUTS:
;
; RETURNS:
;  lls -- LLS structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  QSOS=  -- Structure of statistical QSOs in the sample
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
;   16-Feb-2010 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_llsinit, init


  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'sdss_llsinit, lls, [v1.1]'
      return
  endif 
  
  ;; Initialize
  init = { $
         llsfil: getenv('LLSPAP')+'/SDSS/DR7/Analysis/lls_dr7_stat_LLS.fits', $
         qsofil: getenv('LLSPAP')+'/SDSS/DR7/Analysis/lls_dr7_qsos_sn2050.fits', $
         xzfil: getenv('SDSSPATH')+'/DR7_QSO/xz_val.fits', $
         lX_SLLS: 0.20, $  ;; l(X) for SLLS  (O'Meara et al. 2007, adjusted down)
         vprox: 3000., $
         lambda_LL: 911.7641d, $  ;; Ang
         zmin: 3.3, $  ;; Minimum redshift of the survey
         zem_min: 3.6, $  ;; Minimum redshift of the survey
         zabs_max: 4.4, $  ;; Maximium redshift for LLS analysis
         maxoff: 0.4, $  ;; Maximum redshift offset from zem
         H0: 72., $
         pbins: [ [3.4,3.6], [3.6,3.9], [3.9, 4.2], [4.2,4.4] ], $ ; Proximates
         bins: [ [3.5, 3.65], [3.65,3.9], $
                 [3.9, 4.1], [4.1, 4.4]], $
         all_bins: [ [3.3,3.4], [3.4,3.5], [3.5, 3.65], [3.65,3.9], $
                     [3.9, 4.1], [4.1, 4.4]] $ ;,[4.4,5.0]] $
;               [3.7, 3.85], [3.85, 4.0], [4.0, 4.4], [4.4,5.0] ] $
         }

  return
end
  
