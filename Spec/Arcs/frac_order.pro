;+ 
; NAME:
; frac_order
;     Version 1.1
;
; PURPOSE:
;  For an input array and the order structure defining the echelle
;  footprint, this routine calculates the position of a give pixel
;  along the slit.
;
; CALLING SEQUENCE:
;  frac = frac_order( ordr_str, xstart, ywave, ocen=, ncoeff=)
;
; INPUTS:
;  ordr_str -- Order structure which describes the echelle footprint
;  xstart   -- x values of the pixels
;  ywave    -- Wavelength values of the pixels
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
;   x_slitflat, mike, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   31-Jul-2005 Written by SB
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function frac_order, ordr_str, xstart, ywave, ocen=ocen, ncoeff=ncoeff

      if NOT keyword_set(ncoeff) then ncoeff=5L
      
      nrow = n_elements(ordr_str.lhedg)
      ycol = dindgen(nrow)
      oo = 0.5*(ordr_str.lhedg + ordr_str.rhedg)

      yleft = x_qckwav(ordr_str.lhedg - oo, ycol, ordr_str.arc_m)
      yright = x_qckwav(ordr_str.rhedg - oo, ycol, ordr_str.arc_m)
      all_left  = poly(ywave, poly_fit(yleft, ordr_str.lhedg, ncoeff))
      all_right = poly(ywave, poly_fit(yright, ordr_str.rhedg, ncoeff))

      if NOT keyword_set(ocen) then ocen = 0.5 * (all_right + all_left)

      slit_frac   = 2.0* (xstart - ocen) / (all_right - all_left)

      return, slit_frac
end


