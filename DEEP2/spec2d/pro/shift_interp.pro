function shift_interp,spec,move,spline=spline
;+
; NAME:
;       SHIFT_INTERP
;
; PURPOSE:
;       Shift an array by a non-integer amount 
;       using interpolation.
;
; CALLING SEQUENCE:
;
;       RESULT = SHIFT_INTERP( SPEC, MOVE, /SPLINE )
;
; INPUTS:
;
;       SPEC:   Array to be shifted
;       MOVE:   Amount to shift. If SPEC is 1-D, then MOVE is a
;       scalar. If SPEC is 2-D and if MOVE is a scalar, shift is applied
;       only in X-direction. If MOVE is a 2-D array, then shift is
;       applied in both directions.
;
; OUTPUTS:
;
;       RESULT:  Shifted array
;
; KEYWORD PARAMETERS:
;
;       SPLINE:  Use spline rather than linear interpolation. For 1-D
;       only.
;
; SIDE EFFECTS:
;
;       Don't use this with excessively noisy
;       data lest you interpolate and therefor exacerbate
;       the noise.
;
; EXAMPLE:
;  
; Create a gaussian and shift it to the right 
; by 0.5
;
;     IDL> x = findgen(81)-40
;     IDL> g = exp(-x^2/10.^2)
;     IDL> plot,x,g,/xs
;     IDL> gs = shift_interp(g,0.5)
;     IDL> oplot,x,gs,lines=3
;     IDL> gsreal = exp(-(x-0.5)^2/10.^2)
;     IDL> print,stdev(gsreal[1:*] - gs[1:*])/max(gs)
;            1.0892190e-06
;
; MODIFICATION HISTORY:
; 02.08.2003   Written by JohnJohn
; 04.28.2003   JJ - Fixed problem where nothing would be returned for
; integer shifts. Now handles integer shifts as if it were SHIFT.PRO
; 05.02.2003   JJ - Now uses faster spline method. SPLINE.PRO is
; slooowwwww...
; 06.04.2003   JJ - Now uses linear interpolation by default and spline
; interpolation as a keyword option.
; 06.09.2003   JJ - Corrected mistake by repacing IF bigmove GT 0 with
;                   IF abs(bigmove) GT 0
; 11.02.2005   JJ - Added 2-D shifting for my imaging peeps.
; 06.14.2010   JJ - Fixed 2-D shifting bug to apply the fractional offset in correct direction. Thanks to Jared Rand for finding the bug.
;
;-
on_error,2                      ;If broke, return to sender.

if total(abs(move)) eq 0 then return,spec
fracmove = (move mod 1)
sz = size(spec)
if sz[0] gt 2 then message,'Input array can have no more than 2 dimensions',/io
if sz[0] eq 1 then begin
    if fracmove[0] ne 0 then begin ; do fractional part of shift
        if keyword_set(spline) then begin ;spline onto fractional scale
            x = indgen(n_elements(spec))
            xnew = x + fracmove
        ;;; spl_interp is much faster than spline.
            specnew = spl_interp(xnew,spec,spl_init(xnew,spec),x,/double)
        endif else if fracmove gt 0 then $ ;linear interpoltation
          specnew = (1-fracmove)*spec + fracmove*shift(spec,1) else $
          specnew = (1-abs(fracmove))*spec + abs(fracmove)*shift(spec,-1)
    endif else specnew = spec
endif else begin
    ;;; 2-D fractional shifting
   x = indgen(sz[1]) - fracmove[0]
   if n_elements(move) gt 1 then $
     y = indgen(sz[2]) - fracmove[1] $
   else y = indgen(sz[2])
   specnew = interpolate(spec, x, y, /gr)
endelse
bigmove = fix(move)             ;do integer part of shift 
if total(abs(bigmove)) gt 0 then $
  specbm = shift(specnew, bigmove) $
else specbm = specnew

return,specbm
end
