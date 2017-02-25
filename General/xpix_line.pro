;+ 
; NAME:
; xpix_line
;   Version 1.1
;
; PURPOSE:
;  Finds all pixels within a ray (or line) of specified width.
;
; CALLING SEQUENCE:
;   pval = xpix_line(x0, y0, x1, y1, wid, xmax, ymax, /NOZERO, COUNT=)
;
; INPUTS:
;   x0 -- x position of one point on the ray
;   y0 -- y position of one point on the ray
;   x1 -- x position of one point on the ray
;   y1 -- y position of one point on the ray
;   wid -- Width of the line
;
; RETURNS:
;  pval -- pixels in the circle:  array[2,npt]
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOZERO -- Disallow value of 0 (minimum is 1 instead)
;  COUNT=  -- Number of the pixels in the ray
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function xpix_line, x0, y0, x1, y1, wid, xmax, ymax,  NOZERO=nozero, $
                    COUNT=count

; Finds all pixels (integer values >= 0) within a circle

; 
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
      'arrpix = xpix_line(x0, y0, x1, y1, width, xmax, ymax, /NOZERO, COUNT=)'
    return, -1
  endif 

  if not keyword_set( NOZERO ) then nozero = 0 else nozero = 1

  xmin = nozero 
  ymin = nozero 

; b, m

  m = float(y1-y0)/float(x1-x0)
  b = float(y0) - m*float(x0)

; Check on endpoints

  if( abs(x0-xmin) LT 10 ) then begin
     if( abs(x1-xmax) GT 10 AND abs(y1-ymax) GT 10 AND abs(y1-ymin) GT 10) $
       then begin
         xmax = x1
         if(m LT 0) then ymin = y1 else ymax = y1
     endif
  endif else begin 
      if( abs(y0-ymax) GT 10 AND abs(y0-ymin) GT 10 ) then begin
          xmin = x0
          if (m GT 0) then ymin = y0 else ymax = y0
      endif 
      ; Other endpoint
      if( abs(x1-xmax) GT 10 AND abs(y1-ymax) GT 10 AND abs(y1-ymin) GT 10) $
        then begin
          xmax = x1
          if(m LT 0) then ymin = y1 else ymax = y1
      endif
  endelse
      
  totarr = lonarr(2, long(xmax-xmin+1)*long(ymax-ymin+1))


; Key on b

  count = 0L
  if( b GT 0 AND b LT ymax ) then begin
      for i=xmin, xmax do begin
          y = round(m*float(i) + b)
          ym = (y - wid) > ymin
          yx = (y + wid) < ymax
          for j=ym, yx do begin
              totarr[0,count] = i
              totarr[1,count] = j
              count = count + 1
          endfor
      endfor
  endif else begin
;  Error check
      if (m GT 0 AND b GT ymax) then begin
          print, 'Probelm in xpix_line'
          return, -1
      endif

      my = float(x1-x0)/float(y1-y0)
      by = float(x0) - my*float(y0)
      
      for j=ymax,ymin,-1 do begin
          x = round(my * float(j) + by)
          xm = (x - wid) > xmin
          xx = (x + wid) < xmax
          for i=xm,xx do begin
              totarr[0,count] = i
              totarr[1,count] = j
              count = count + 1
          endfor
      endfor
  endelse
          
; Return

  if count GT 0 then return, totarr[0:1,0:count-1] else return, -1

end
          
