;+ 
; NAME:
; xgetx_plt
;    Version 1.1
;
; PURPOSE:
;  Returns the x-value corresponding to the x-pixel on a GUI.
;   Requires pos, xymnx or the appropriate Structure.
;
; CALLING SEQUENCE:
;  xval = xgetx_plt(xtmp, pos, xymnx, size, /STRCT)
;   
; INPUTS:
;  xtmp -- x-pixel value
;  pos  -- Fraction of plot window covered by plot (2 element array)
;  xymnx -- x-y limits of plot window (4 element array: x0,y0,x1,y1)
;
; RETURNS:
;   xval -- 
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
;   uniq = x_uniqstr( lbls, count=count)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function xgetx_plt, xtmp, pos, xymnx, size, STRCT=strct, LOG=log

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x = xgetx_plt(xpix,pos,xymnx,size)'
    print,'or Syntax - ' + $
             'x = xgetx_plt(strc, /STRCT)'
    return, -1
  endif 

; Allow for structure
  if not keyword_set(LOG) then log = 0
  if keyword_set( STRCT ) then begin
      xpix = xtmp.xcurs
      pos = xtmp.pos
      xymnx = xtmp.xymnx
      idx = tag_exist(xtmp, 'winsize', /top_level)
      if idx EQ 0 then size = xtmp.size else size = xtmp.winsize
      if tag_exist(xtmp, 'flg_logx', /top_level) then log = xtmp.flg_logx
  endif else xpix = xtmp
  
; Returns the data value from the pixel value
  if keyword_set(LOG) then begin
      m = (alog10(xymnx[2]) - alog10(xymnx[0])) / (size[0] * (pos[2] - pos[0]))
      b = alog10(xymnx[2]) - m * size[0] * pos[2]
      val = 10^(m*xpix + b)
  endif else begin
      m = (xymnx[2] - xymnx[0]) / (size[0] * (pos[2] - pos[0]))
      b = xymnx[2] - m * size[0] * pos[2]
      val = m*xpix + b
  endelse
  return, val
end
