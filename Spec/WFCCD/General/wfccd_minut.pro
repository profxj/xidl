;+ 
; NAME:
; wfccd_minut   Version 1.0
;
; PURPOSE:
;    Creates Bias frame given structure
;
; CALLING SEQUENCE:
;   
;   wfccd_minut, struct, mask_id, SVOV=svov
;
; INPUTS:
;   struct -- wfccd_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   flat - fits file in the dir Flats named 'Flat_##.fits'
;                 where ## is the mask_id value
;   VAR  - Variance in the flat (in electrons)
;
; OPTIONAL KEYWORDS:
;   SVOV - save ov files
;   NOFITS - No FITS output
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_minut, nght1_strct, mask_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function wfccd_minut, wfccd, scalar, vector, MN=mn

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'imn = wfccd_minut( struct, scalar, vector ) [v1.0]'
      return, -1
  endif 
  
;  Optional Keywords
  
;  Check scalar 
  if n_elements(scalar) GT 1 then begin
      print, 'wfccd_minut: Scalar must be scalar!'
      return, -1
  endif

;  Parse scalar
  hr_s = long(strmid(wfccd[scalar].ut,0,2))
  min_s = long(strmid(wfccd[scalar].ut,3,2))
  sec_s = float(strmid(wfccd[scalar].ut,6,4))
  all_s = hr_s*3600. + min_s*60. + sec_s

;  Vector
  hr = long(strmid(wfccd[vector].ut,0,2))
  min = long(strmid(wfccd[vector].ut,3,2))
  sec = float(strmid(wfccd[vector].ut,6,4))
  all = hr*3600. + min*60. + sec

;  Find min
  diff = abs(all_s-all)
  ; Allow for flipping of the clock at 24hr
  a = where( diff GT 20*3600., na)
  if na NE 0 then diff[a] = abs(diff[a]-24*3600.)

  mn = min(diff, imn)

; Return
  return, imn

end
