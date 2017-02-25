;+ 
; NAME:
;  x_mkwvarray   
;   Version 1.1
;
; PURPOSE:
;    Generate a wavelength array given some input
;
; CALLING SEQUENCE:
;   
;   wave = x_setwave(head, ntot)
;
; INPUTS:
;   del      - Dispersion in km/s
;    w0      - Starting wavelength
;    w1      - Ending wavelength
;
; RETURNS:
;   wave     - Wavelength array
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
;   wave = x_setwave(head, 1000L)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP 
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_mkwvarray, del, w0, w1, DWAVE=dwave, NPIX=npix

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'wave = x_mkwvarray(del, w0, [w1]) [v1.1]'
    return, -1
  endif 

  c = x_constants()
  ;;  Optional Keywords
  if not keyword_set(DWAVE) then begin
      ;; Assumes constant del in km/s
      cdelt = del*1d5 / c.c / alog(10.)    
      npix = alog10(w1 / w0)  / cdelt + 1
      wave = w0 * 10.d^(dindgen(npix)*cdelt)
  endif else stop

  return, wave

end
