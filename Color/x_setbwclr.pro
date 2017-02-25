;+ 
; NAME:
; x_setbwclr
;
; PURPOSE:
;    Creates an array of B&W shades useful for converting index
;    numbers into 'color'
;
; CALLING SEQUENCE:
;   
;   xcolors = x_setclrs()
;
; INPUTS:
;   [nclr]  - Optional input setting the array size (Default is 12)
;
; RETURNS:
;   xcolors - String Array of B&W shades
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
;  xcolors = x_setbwclr()
;
; PROCEDURES CALLED:
;   getcolor (Coyote package)
;
; REVISION HISTORY:
;   25-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function x_setbwclr, nclr, PSFILE=psfile, SCL=scl

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'clr = x_setbwclr(nclr, /PSFILE) [v1.0]'
    return, -1
  endif 
;
  if not keyword_set( SCL ) then  scl = 1.

;  Get one color

  resolve_routine, 'getcolor', /no_recompile, /is_function

  if not keyword_set(PSFILE) then allclr = lonarr(nclr) $
  else allclr = lonarr(nclr,3)

  val = round(255.*findgen(nclr)/(float(nclr)*scl))
  for qq=0L,nclr-1 do begin
      vec = replicate(val[qq],3)
      if not keyword_set( PSFILE ) then allclr[qq] = color24(vec) $
      else allclr[qq,*] = vec
  endfor

  return, allclr

end
