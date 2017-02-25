;+ 
; NAME:
; iraf_satt 
;    Version 1.1
;
; PURPOSE:
;    Create a set of IRAF directives to mask out a satellite
;  trail in a direct image
;
; CALLING SEQUENCE:
;   iraf_satt, x1, y1, x2, y2, mskfil, dx=, dy=, SATFIL=, YMAX=, 
;     XMAX=, XMIN=
;
; INPUTS:
;   x1,y1  -- x,y positions at one point (end) of the trail
;   x2,y2  -- x,y positions at another point (end) of the trail
;   mskfil -- Name of mask file to edit 
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   dx  -- Width of satellite trail in x dimension
;   dy  -- Width of satellite trail in y dimension
;  satfil -- Name of satellite file
;  ymax -- Maximum y value for the mask
;  xmax -- Maximum x value for the mask
;  xmin -- Minimum x value for the mask
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   iraf_satt, 10., 50., 100., 90., dx=5
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Nov-2001 Written by JXP
;   06-Feb-2002 Revised to allow FILL (JXP)
;-
;------------------------------------------------------------------------------
pro iraf_satt, x1,y1,x2,y2, mskfil, dx=dx, dy=dy, SATFIL=satfil, YMAX=ymax,$
               XMAX=xmax, XMIN=xmin

  if  N_params() LT 5  then begin 
      print,'Syntax - ' + $
        'iraf_satt, x1, y1, x2, y2, mskfil, dx=, dy=, satfil=, YMAX=, XMAX='
      print, '       XMIN= [v1.1]'
      return
  endif 

  if not keyword_set(dx) and not keyword_set(dy) then begin
      print, 'Must set dx or dy'
      return
  endif

  if not keyword_set( SATFIL ) then satfil='satt.cl'
  if not keyword_set( TMPFIL ) then tmpfil='tmp_sat.pl'

  ;; Deal with tmpfil
  a = findfile(tmpfil,count=na)
  if na NE 0 then spawn, 'rm -f '+tmpfil
 
  ;; dx
  if keyword_set(dx) then begin
      ;; Get slope, b  x=my+b
      m = float(x2-x1)/float(y2-y1)
      b = x1 - m*y1

      ;; Create imexpr
      close, /all
      openw, 1, 'expr'
      expr = '( (abs('+strtrim(m,2)+'*J +'+strtrim(b,2)+'-I) < '+$
        strtrim(long(dx),2)+')'
      ;; YMAX
      if keyword_set(YMAX) then $
        expr = expr + ' && (J < '+strtrim(long(ymax),2)+')'

      ;; CLOSE EXPR
      expr = expr+') ? 1 : a'
      printf, 1, expr
      close, /all
      
      ;; Create satfil
      openw, 1, satfil
      printf, 1, 'imdel '+tmpfil
      printf, 1, 'imexpr @expr '+tmpfil+' a='+mskfil
      printf, 1, 'imdel '+mskfil
      printf, 1, 'imrename '+tmpfil+' '+mskfil
      close, /all
  endif

  ;; dy
  if keyword_set(dy) then begin
      ;; Get slope, b  y=mx+b
      m = float(y2-y1)/float(x2-x1)
      b = y1 - m*x1

      ;; Create imexpr
      close, /all
      openw, 1, 'expr'
      expr = '( (abs('+strtrim(m,2)+'*I +'+strtrim(b,2)+'-J) < '+$
        strtrim(long(dy),2)+')'
      ;; YMAX
      if keyword_set(YMAX) then $
        expr = expr + ' && (J < '+strtrim(long(ymax),2)+')'
      ;; XMAX
      if keyword_set(XMAX) then $
        expr = expr + ' && (I < '+strtrim(long(xmax),2)+')'
      ;; XMIN
      if keyword_set(XMIN) then $
        expr = expr + ' && (I > '+strtrim(long(xmin),2)+')'

      ;; CLOSE EXPR
      expr = expr+') ? 1 : a'
      printf, 1, expr
      close, /all
      
      ;; Create satfil
      openw, 1, satfil
      printf, 1, 'imdel '+tmpfil
      printf, 1, 'imexpr @expr '+tmpfil+' a='+mskfil
      printf, 1, 'imdel '+mskfil
      printf, 1, 'imrename '+tmpfil+' '+mskfil
      close, /all
  endif

  return
end

      

      
