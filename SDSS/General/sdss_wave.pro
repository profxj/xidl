;+ 
; NAME:
; sdss_wave
;    Version 1.1
;
; PURPOSE:
;    Create a wavelength array for SDSS data
;
; CALLING SEQUENCE:
;  sdss_wave, wave, pxiels=, iwave=
;
; INPUTS:
;
; RETURNS:
;  wave -- Wavelength array
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
;   Written by SHF
;-
;------------------------------------------------------------------------------
pro sdss_wave, wave, pixels=pixels, iwave=iwave
if not keyword_set(pixels) then pixels= 3850L
if not keyword_set(iwave) then iwave= 3.5793d
wave= 10^(iwave + findgen(pixels)*(0.0001))-10^(iwave)+10^(iwave)
end
