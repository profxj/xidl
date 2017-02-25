;+ 
; NAME:
; x_fndreg
;   Version 1.1
;
; PURPOSE:
;    Finds data points within a number of regions.  
;    [This is a rather difficult routine to use]
;
; CALLING SEQUENCE:
;   
;   pts = x_fndreg(xdat, reg)
;
; INPUTS:
;   xdat - Data
;   reg  - Regions:   fltarr(N,2)
;
; RETURNS:
;   pts  - Points in the regions
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   NPNT -  Number of points within the regions
; COMMENTS:
;
; EXAMPLES:
;   pts = x_fndreg( findgen(1000), [15., 25.])
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   21-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------

function x_fndreg, xdat, reg, NPNT=npnt

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'pts = x_fndreg(xdat, reg, NPNT=)  [V1.1]'
    return, -1
  endif 

; Check to make sure no regions overlap

  sz = size(reg, /dimensions)
  nreg = sz[0]
  dumi = lindgen(nreg)
  
  if nreg NE 1 then begin
      for i=0L, nreg-1 do begin
          noti = where(dumi NE i)
          bla = where(reg[noti,*] LT reg[i,1] AND reg[noti,*] GT reg[i,0], count)
          if count NE 0 then begin
              print, 'Regions may not overlap!'
              npnt = 0
              return, -1
          endif
      endfor
  endif

; Find the points

  flg_reg = 0
  pts = -1
  for i=0L,nreg-1 do begin
      dum = where(xdat GE reg[i,0] AND xdat LE reg[i,1], count)
      if count NE 0 then begin
          if flg_reg EQ 0 then begin
              pts = dum 
              flg_reg = 1
          endif else pts = [pts, dum]
      endif
  endfor

  if pts[0] NE -1 then npnt = n_elements(pts) else npnt = 0

  return, pts

end

