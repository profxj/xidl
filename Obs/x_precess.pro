;+ 
; NAME:
; x_fndchrt   
;    Version 1.1
;
; PURPOSE:
;  Precess RA and DEC
;
; CALLING SEQUENCE:
;  
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   imsize - Arcmin of image (default is 5')
;
; OPTIONAL OUTPUTS:
;  OUTDIR=  -- Name of output directory
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  precess
;
; REVISION HISTORY:
;   21-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

function x_precess, radec, oldepoch, newepoch

  if  N_params() LT 3 then begin 
      print,'Syntax - ' + $
        'newrad = x_precess(radec, oldepoch, newepoch) [v1.0]'
      return, -1
  endif 

  x_radec, radec[0], radec[1], rad, decd

  precess, rad, decd, oldepoch, newepoch
  ;; For labeling of J2000
  x_radec, ras, decs, rad, decd, /flip

  return, [ras, decs]
end


