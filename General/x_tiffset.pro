;+ 
; NAME:
; x_tiffset
;   Version 1.1
;
; PURPOSE:
;  Sets (or unsets) plotting variables for nice looking screen plots 
;  to then be used in laptop presentations.
;
; CALLING SEQUENCE:
;   
;   x_tiffset, /UNSET
;
; INPUTS:
;   
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /UNSET -- Unset the plotting values
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Jun-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_tiffset, UNSET=unset

;
;  Optional Keywords
  if not keyword_set( UNSET ) then begin
      set_plot, 'x'
      !p.font = 1
      device, set_font='Helvetica Bold', /tt_font
;      device, set_font='Courier', /tt_font
;      device, set_font='Courier Bold', /tt_font
;      device, set_font='Times Bold', /tt_font
      !p.thick = 3
      !x.thick = 3
      !y.thick = 3
      !p.charthick = 3
;      device, /times,isolatin=1
  endif else begin
      !p.thick = 1
      !p.charthick = 1
      !p.font = -1
      !x.thick = 1
      !y.thick = 1
  endelse

end
