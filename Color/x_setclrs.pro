;+ 
; NAME:
; x_setclrs
;
; PURPOSE:
;    Creates an array of colors useful for loops of undetermined
;    length
;
; CALLING SEQUENCE:
;   
;   xcolors = x_setclrs, [nclr]
;
; INPUTS:
;   [nclr]  - Optional input setting the array size (Default is 10)
;
; RETURNS:
;   xcolors - String Array of colors
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
;  xcolors = x_setclrs
;
; PROCEDURES CALLED:
;   getcolor (coyote)
;
; REVISION HISTORY:
;   25-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function x_setclrs, nclr, WHITE=white

;  if (N_params() LT 3) then begin 
;    print,'Syntax - ' + $
;             'getabnd, nam, Zval, abnd, [flag=]'
;    return
;  endif 

;
  if not keyword_set( nclr ) then    nclr = 12

  if nclr GT 12 then message, 'x_setclrs: Not set for nclr > 10'

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
          3: allclr[q] = clr.green
          4: allclr[q] = clr.orange
          5: allclr[q] = clr.purple
          6: allclr[q] = clr.cyan
          7: allclr[q] = clr.gray
          8: allclr[q] = clr.brown
          9: allclr[q] = clr.yellow
          10: allclr[q] = clr.skyblue
          11: allclr[q] = clr.darkgreen
          else : stop
      endcase
  endfor

  return, allclr

end
