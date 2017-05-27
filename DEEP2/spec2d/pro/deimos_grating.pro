
;+
; NAME:
;   deimos_grating
;
; PURPOSE:
;   compute grating angle and rule from FITS header
;
; CALLING SEQUENCE:
;   deimos_grating, header, g_rule, grangle, lambda_c
; 
; INPUTS:
;   header  - FITS header
;
; OUTPUTS:
;   g_rule  - rule spacing [lines/mm]
;   grangle - rule spacing [degrees]
;   lambda_c - central wavelength [angstroms]
;
; MODIFICATION HISTORY:
;   2002-Jun-04  - DPF, MD
;   2002-Jul-10  - doesn't crash on imaging files - DPF
;-

pro deimos_grating, header, g_rule, grangle, lambda_c

  grating_number = fxpar(header, 'GRATEPOS')
  grating_name = fxpar(header, 'GRATENAM' )
  
  if strmid(grating_name, 0, 6) eq 'Mirror' then g_rule = 0 else $
    g_rule = double(grating_name)

  if abs(g_rule-1200.d0) lt 0.5 then g_rule=1200.06d0
  if abs(g_rule-831.d0) lt 2 then g_rule=831.90d0

  
  if grating_number eq 3 then begin
     rawpos = fxpar(header, 'G3TLTRAW')
     grangle = (rawpos+29094)/2500. ;this is subject to change!!
; shouldn't the offset be a FITS keyword??
     lambda_c = fxpar(header, 'G3TLTWAV')
     
  endif else begin
     rawpos = fxpar(header, 'G4TLTRAW')
     grangle = (rawpos+40934)/2500. ;subject to change!!
     lambda_c = fxpar(header, 'G4TLTWAV')
  endelse

  return
end

