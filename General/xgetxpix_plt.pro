;+ 
; NAME:
; xgetx_plt
;    Version 1.1
;
; PURPOSE:
;  Returns the x-pixel value corresponding to the x-value.
;   Requires pos, xymnx or the appropriate Structure.
;
; CALLING SEQUENCE:
;  xpix = xgetx_plt(xtmp, pos, xymnx, size, /STRCT)
;   
; INPUTS:
;  xtmp -- x-pixel value
;  pos  -- Fraction of plot window covered by plot (2 element array)
;  xymnx -- x-y limits of plot window (4 element array: x0,y0,x1,y1)
;
; RETURNS:
;   xpix -- 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STRCT -- xtmp contains a structure with the relevant tags
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
;   17-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

; Returns the pixel value from the data value
function xgetxpix_plt, x, pos, xymnx, size

  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'xpix = xgetxpix_plt(x,pos,xymnx,size)'
    return, -1
  endif 


  m = size[0] * (pos[2] - pos[0]) / (xymnx[2] - xymnx[0]) 
  b = size[0] * pos[2] - m * xymnx[2]
  return, m*x + b
end

