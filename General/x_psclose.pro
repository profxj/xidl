;+ 
; NAME:
; x_psclose
;   Version 1.1
;
; PURPOSE:
;    Call ps_close and reset a number of default values.  Point window
;   at X-windows.
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
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  ps_close
;  set_plot
;
; REVISION HISTORY:
;   09-Feb-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_psclose

;  Optional Keywords


  ps_close, /noprint, /noid
  set_plot, 'x'
  device, decompose=1
  !p.thick = 1
  !p.charthick = 1
  !p.font = -1
  !x.thick = 1
  !y.thick = 1

end
