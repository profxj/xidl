;+ 
; NAME:
; xgety_plt
;    Version 1.1
;
; PURPOSE:
;  Returns the y-pixel corresponding to the y-value on a GUI.
;   Requires pos, xymnx or the appropriate Structure.
;
; CALLING SEQUENCE:
;  yval = xgety_plt(ytmp, pos, xymnx, size, /STRCT)
;   
; INPUTS:
;  ytmp -- y-pixel value
;  pos  -- Fraction of plot window covered by plot (2 element array)
;  xymnx -- x-y limits of plot window (4 element array: x0,y0,x1,y1)
;
; RETURNS:
;   yval -- 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STRCT -- ytmp contains a structure with the relevant tags
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   uniq = x_uniqstr( lbls, count=count)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function xgetypix_plt, y, pos, xymnx, size

  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'ypix = xgetypix_plt(y,pos,xymnx,size)'
    return, -1
  endif 

; Returns the data value from the pixel value

  m = size[1] * (pos[3] - pos[1]) / (xymnx[3] - xymnx[1]) 
  b = size[1] * pos[3] - m * xymnx[3]
  return, m*y + b
end
