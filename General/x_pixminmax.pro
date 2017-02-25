;+ 
; NAME:
; x_pixminmax
;  Version 1.1
;
; PURPOSE:
;  Find pixels in the wavelength array corresponding to the redshift
;  and rest wavelength.
;
; CALLING SEQUENCE:
;   
;   x_pixminmax, wave, wrest, zabs, [vmin, vmax], PIXMIN=, PIXMAX=, VELO=
;
; INPUTS:
;  wave -- Wavlength array
;  wrest -- Rest wavelength of transition
;  zabs -- Absorption redshift
;  [vmin,vmax] -- Velocity minimum and maximum
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  PIXMIN= -- Pixel value corresponding to VMIN
;  PIXMAX= -- Pixel value corresponding to VMAX
;  VELO= -- Velocity array
;
; COMMENTS:
;
; EXAMPLES:
;   x_pixminmax, wave, 1808.0126d, 1.5, -100., 50., PIXMIN=pixmn
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
    print, '              VELO= [v1.1]'
    return
  endif 

  spl=2.9979e5

;; Create VELO

  velo = (wave-wrest*(1+zabs))*spl/( wrest*(1+zabs) )
  ;; PIXMIN, PIXMAX
  if arg_present(PIXMIN) AND (n_elements(VMIN) NE 0) then $
    mn = min(abs(velo-vmin),pixmin)
  if arg_present(PIXMAX) AND (n_elements(VMAX) NE 0) then $
    mn = min(abs(velo-vmax),pixmax)


  return
end

