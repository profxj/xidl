;+ 
; NAME:
;  x_tvx
;   Version 1.1
;
; PURPOSE:
;    Returns the x value of the image from the pixel value
;
; CALLING SEQUENCE:
;   
;   x = x_tvx(tv, /INTG)
;
; INPUTS:
;   tv -- TV structure  (tvstrct)
;
; RETURNS:
;   x  -- x value from the image from the pixel value
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /INTG  -- Returns an integer value
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; 
; PROCEDURES/FUNCTIONS CALLED:

; REVISION HISTORY:
;   
;
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_tvx, tv, INTG=intg

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x = x_tvx(tv, /INTG) '
    return, -1
  endif 

; Returns the data value from the pixel value

  m = (tv.xymnx[2] - tv.xymnx[0]) / (tv.winsize[0] * (tv.pos[2] - tv.pos[0]))
  b = tv.xymnx[2] - m * tv.winsize[0] * tv.pos[2]

  if keyword_set( INTG ) then return, round(m*tv.xcurs + b) $
    else return, m*tv.xcurs + b 

end
