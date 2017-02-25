;+ 
; NAME:
; x_psopen
;   Version 1.1
;
; PURPOSE:
;  Calls ps_open and sets a number of default plotting values.  
;
; CALLING SEQUENCE:
;   
;   x_psopen, fil, /MAXS, /PORTRAIT
;
; INPUTS:
;   img       - Fits file or data
;
; RETURNS:
;   dat       - Data in fits file or the input data
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  FSCALE      - Data is float
;
; OPTIONAL OUTPUTS:
;  HEAD        - Header
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Feb-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_psopen, psfile, MAXS=maxs, _EXTRA=extra

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_psopen, psfile, /MAXS, _EXTRA= [V1.1]'
    return
  endif 

;  Optional Keywords
  set_plot, 'x'
  !p.font = 0
  device, decompose=0
  ps_open, filename=PSFILE, /color, bpp=8, MAXS=maxs, _EXTRA=extra
  !p.thick = 6
  !x.thick = 6
  !y.thick = 6
  !p.charthick = 3
  device, /times,isolatin=1

end
