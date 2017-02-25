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
function mtham_trans, wave
;
;  Extintion at Mauna Kea, taken from CFHT Bulletin, 19, 16 (1988). 
;  I/O - communicate extinction at this wavelenth
;..............................................................................
;
  ;; Extinction file
  longslit_dir = getenv('LONGSLIT_DIR')
  extinctfile = longslit_dir + '/calib/extinction/mthamextinct.dat'
  readcol, extinctfile, wave_ext, mag_ext, format = 'F,F' , /silent

  ;; Interpolate
  wave_min = min(wave_ext, jmin)
  wave_max = max(wave_ext, jmax)
  mag_ext1 = dblarr(n_elements(wave))
  ext = dblarr(n_elements(wave)) + 1.0D
  inds = WHERE(wave GE wave_min AND wave LE wave_max)
  mag_ext1[inds] = interpol(mag_ext, wave_ext, wave[inds])
  linds = WHERE(wave LT wave_min, nl)
  IF nl GT 0 THEN mag_ext1[linds] = mag_ext[jmin]
  rinds = WHERE(wave GT wave_max, nr)
  IF nr GT 0 THEN mag_ext1[rinds] = mag_ext[jmax]

  return, mag_ext1

end

