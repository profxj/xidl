;+ 
; NAME:
; xpix_circ
;   Version 1.1
;
; PURPOSE:
;  Finds all pixels within a circle on an image.
;
; CALLING SEQUENCE:
;   pval = xpix_circ(x0, y0, radius, /NOZERO, MAXX=, MAXY=, COUNT=)
;
; INPUTS:
;   x0 -- Circle x position
;   y0 -- Circle y position
;   radius -- Radius of the circle
;
; RETURNS:
;  pval -- pixels in the circle:  array[2,npt]
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOZERO -- Disallow value of 0 (minimum is 1 instead)
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
;   17-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function xpix_circ, x0, y0, radius, NOZERO=nozero, MAXX=maxx, MAXY=maxy, $
                    COUNT=count

; Finds all pixels (integer values >= 0) within a circle

; 
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'arrpix = xpix_circ(x0, y0, radius, /NOZERO, MAXX=, MAXY=, COUNT=)'
    return, -1
  endif 

  if not keyword_set( NOZERO ) then nozero = 0 else nozero = 1

  xmin = nozero > round(x0 - radius) 
  if keyword_set( MAXX ) then xmax = maxx < round(x0 + radius) $
    else xmax = round(x0 + radius)
  ymin = nozero > round(y0 - radius) 
  if keyword_set( MAXY ) then ymax = maxy < round(y0 + radius) $
    else ymax = round(y0 + radius)

; Return if the circle isnt on the CCD

  if (xmax LT xmin OR ymax LT ymin) then begin
      count = 0
      return, -1
  endif

; Loop

  radsq = radius^2
  totarr = intarr(2, long(xmax-xmin+1)*long(ymax-ymin+1))
  count = 0L
  for i=xmin, xmax do begin
      for j=ymin,ymax do begin
          if( (i-x0)^2 + (j-y0)^2 LT radsq) then begin
              totarr[0,count] = i
              totarr[1,count] = j
              count = count + 1
          endif
      endfor
  endfor

; Return

  if count GT 0 then return, totarr[0:1,0:count-1] else return, -1

end
          
