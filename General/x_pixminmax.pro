;+ 
; NAME:
; x_pixminmax
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;   x_pixminmax, ra, dec, rad, decd, /ARCS, /FLIP
;
; INPUTS:
;   ra, dec    - RA and DEC in in RR:RR:RR.R -DD:DD:DD.D format 
;                 Colons are required as separators
;   rad, decd  - RA and DEC in decimal degrees (double)
;
; RETURNS:
;
; OUTPUTS:
;   ra, dec    - RA and DEC in in RR:RR:RR.R -DD:DD:DD.D format 
;                 Colons are required as separators
;   rad, decd  - RA and DEC in decimal degrees (double)
;
; OPTIONAL KEYWORDS:
;   ARCS - Outputs in arcseconds
;   FLIP - Gives RA and DEC from decimal RA,DEC
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_pixminmax, '21:12:23.1', '-13:13:22.2', rad, decd
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Jun-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_pixminmax, wave, wrest, zabs, vmin, vmax, $
                 PIXMIN=pixmin, PIXMAX=pixmax, VELO=velo

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'x_pixminmax, wave, wrest, zabs, [vmin,vmax], PIXMIN=, PIXMAX=, '
    print, '              VELO= [v1.0]'
    return
  endif 

  spl=2.9979e5

;; Create VELO

  velo = (wave-wrest*(1+zabs))*spl/( wrest*(1+zabs) )
  ;; PIXMIN, PIXMAX
  if arg_present(PIXMIN) AND keyword_set(VMIN) then mn = min(abs(velo-vmin),pixmin)
  if arg_present(PIXMAX) AND keyword_set(VMAX) then mn = min(abs(velo-vmax),pixmax)


  return
end

