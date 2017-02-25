;+ 
; NAME:
; manuakea_trans
;    Version 1.1
;
; PURPOSE:
;    Extinction of the atmosphere at Mauna Kea
;
; CALLING SEQUENCE:
;  extinct = maunakea_trans(wave)
;
; INPUTS:
;  wave=  -- Wavelengths of interest
;
; RETURNS:
;  extinct  -- Extinction in mag 
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
;  showfits
;  querydss
;
; REVISION HISTORY:
;   27-Oct-2005 Written by JXP based on HIRES S2N code
;-
;------------------------------------------------------------------------------
function maunakea_trans, wave
;
;  Extintion at Mauna Kea, taken from CFHT Bulletin, 19, 16 (1988). 
;  I/O - communicate extinction at this wavelenth
;..............................................................................
;
 xwave = [3000,3100,3200, 3300, 3400., $
          3500,3600,3700, 3800, 3900, $
          4000,4250,4500, 4750, 5000, $
          5250,5500,5750, 6000, 6500, $
          7000,8000,9000,10000,12000 ]
;
 xthru = [4.90, 1.37, 0.82, 0.57, 0.51, $
          0.42, 0.37, 0.33, 0.30, 0.27, $
          0.25, 0.21, 0.17, 0.14, 0.13, $
          0.12, 0.12, 0.12, 0.11, 0.11,$
          0.10, 0.07, 0.05, 0.04, 0.03]

  return, interpol(xthru, xwave, wave)

end

