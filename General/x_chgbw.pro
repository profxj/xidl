;+ 
; NAME:
; x_chgbw
;   Version 1.1
;
; PURPOSE:
;   Sets colors to black and white 
;
; CALLING SEQUENCE:
;   
;   x_chgbw
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
;   x_chgbw
;
; PROCEDURES CALLED:
; getcolor
;
; REVISION HISTORY:
;   19-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------
pro x_chgbw


  ; Color
  clr = getcolor(/load)
  !p.color = clr.black
  !p.backgournd = clr.white

  return
end
