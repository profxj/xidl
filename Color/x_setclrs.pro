;+ 
; NAME:
; x_setclrs
;
; PURPOSE:
;    Creates an array of colors useful for converting index
;    numbers into color
;
; CALLING SEQUENCE:
;   
;   xcolors = x_setclrs()
;
; INPUTS:
;   [nclr]  - Optional input setting the array size (Default is 12)
;
; RETURNS:
;   xcolors - String Array of colors
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /WHITE  - Assumes a black background 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  xcolors = x_setclrs()
;
; PROCEDURES CALLED:
;   getcolor (Coyote package)
;
; REVISION HISTORY:
;   25-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function x_setclrs, nclr, WHITE=white, DARK=dark

;
  if not keyword_set( nclr ) then    nclr = 20
  if not keyword_set( init ) then    init = 0L

  if nclr GT 20 then message, 'x_setclrs: Not set for nclr > 20'

;  Get one color

  clr = getcolor(/load)
  a = clr.yellow
  allclr = replicate(a, nclr)

  for q=0,nclr-1 do begin
      case q of 
          0: begin
              if keyword_set(WHITE) then allclr[q] = clr.white else $
                allclr[q] = clr.black
          end
          1: allclr[q] = clr.blue
          2: allclr[q] = clr.red
          3: begin
             if not keyword_set(DARK) then allclr[q] = clr.green $
              else allclr[q] = clr.darkgreen 
          end
          4: allclr[q] = clr.orange
          5: allclr[q] = clr.purple
          6: allclr[q] = clr.cyan
          7: allclr[q] = clr.darkgray
          8: allclr[q] = clr.darkgreen
          9: allclr[q] = clr.yellow
          10: allclr[q] = clr.skyblue
          11: allclr[q] = clr.brown
          12: allclr[q] = clr.pink
          13: allclr[q] = clr.maroon
          14: allclr[q] = clr.tan
          15: allclr[q] = clr.lightgray
          16: allclr[q] = clr.firebrick
          17: allclr[q] = clr.gray
          18: allclr[q] = clr.cornsilk
          19: allclr[q] = clr.seagreen
          20: allclr[q] = clr.orangered
          else : stop
      endcase
  endfor

  return, allclr

end
