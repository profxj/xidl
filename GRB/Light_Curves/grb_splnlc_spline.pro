;+ 
; NAME:
;  grb_splnlc
;   Version 1.1
;
; PURPOSE:
;    Given a grb structure (magnitude at a time t,  z, spectral slope)
;    calculates the luminsoity at any time 
;
; CALLING SEQUENCE:
;   
;   lum_nu = grb_calclum(grb, nu, [t])
;
; INPUTS:
;     grb -- GRB structure
;     nu  -- Frequency (GRB frame)
;     [t] -- Time (observer frame). Default = t0 in GRB structure
;
; RETURNS:
;   lum_nu= 
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
;   17-Feb-2006 Written by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function grb_splnlc_spline, t

common grb_splinelc, grbsp_nu, grbsp_mag, grbsp_t, grbsp_splin

  val = spl_interp(grbsp_t, grbsp_mag, grbsp_splin, t, /double)
  return, val

end
