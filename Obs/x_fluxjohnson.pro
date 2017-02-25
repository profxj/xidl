;+ 
; NAME:
; x_fluxjohnson   
;    Version 1.1
;
; PURPOSE:
;    Calculates conversion factor for Johnson magnitudes to physical
;    flux by interpolation 
;
; CALLING SEQUENCE:
;  flux_conv = x_fluxjohnson(wave)
;
; INPUTS:
;  wave= -- Wavelength array
;
; RETURNS:
;  flux_conv -- Conversion factor from mag to flux
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
;   30-Aug-2005 Written by JXP based on HIRES S2N code
;-
;------------------------------------------------------------------------------
function x_fluxjohnson, wave

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'flux = x_fluxjohnson(wave) (v1.0)'
      return, -1
  endif 

  xwave = [3062.,3312,3562,3812,4062,4212,4462,4712,4962,5425,5925, $
           6425,6925,7750,8350,11050]
  xflux = [569.,586,590,1326,1770,1707,1530,1356,1257,1054,886,749, $
           641,502,435,263]

  return, interpol(xflux, xwave, wave)

end


