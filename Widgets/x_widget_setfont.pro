;+ 
; NAME:
; x_widget_setfont
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
;  x_screen -- Screen size in pixels
;
; RETURNS:
;   strct --  A structure of fonts 
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

function x_widget_setfont, x_screen

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x = x_widget_slider, strct, [indx], /INIT, /REINIT, /DRAG [v1.0]'
    return, -1
  endif 

  fonts = {widgetfontstrct}
  fonts.scr_xpix = x_screen

  case !d.name of
      'X': begin  ; X windows
          if abs(x_screen - 1400) LT 200 then begin ;; xpix = 1200 to 1600
              fonts.small_font = '7x13'
              fonts.big_font = '9x15'
          endif else if abs(x_screen - 1800) LT 200 then begin ;; 1600 to 2000
              fonts.small_font = '8x13'
              fonts.big_font = '10x20'
          endif else begin
              print, 'Your screen has an x pixel range that is untested'
              print, 'Hopefully the fonts will work ok for you'
              print, 'Continue on...'
              stop
          endelse
      end
      else: stop
  endcase

  return, fonts
end
