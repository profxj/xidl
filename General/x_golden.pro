;+ 
; NAME:
; x_golden
;   Version 1.1
;
; PURPOSE:
;    Uses the NR routine 'golden' to find the x-value where
;    the function has a minimum
;
; CALLING SEQUENCE:
;   
;   min = x_golden(func,a,b,c TOL=)
;
; INPUTS:
;   func - String name of the IDL function
;   a,c  - Values bracketing the minimum
;   b    - Best guess at minimum
;
; RETURNS:
;   min - x-value of function minimum
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   TOL - Fractional Tolerance
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   num = x_golden('func', 0.0, 1.0)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   25-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_golden, func, a, b, c, TOL=tol 

;
  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'line = x_golden(func, a, b, c, TOL=) [v1.1]'
    return, -1
  endif 

;  Optional Keywords
  if not keyword_set( TOL ) then tol = 1.d-5


; Giddeup
  r_gld = 0.61803399d
  c_gld = 1.d - r_gld

  x0 = a
  x3 = c

  if abs(c-b) GT abs(b-a) then begin
      x1 = b
      x2 = b+c_gld*(c-b)
  endif else begin
      x2 = b
      x1 = b-c_gld*(b-a)
  endelse

  f1 = call_function(func, x1)
  f2 = call_function(func, x2)

  while abs(x3-x0) GT TOL*(abs(x1)+abs(x2)) do begin
      if f2 LT f1 then begin
          x0 = x1
          x1 = x2
          x2 = r_gld*x1 + c_gld*x3
          f1 = f2
          f2 = call_function(func, x2)
      endif else begin
          x3 = x2
          x2 = x1
          x1 = r_gld*x2 + c_gld*x0
          f2 = f1
          f1 = call_function(func, x1)
      endelse
  endwhile

  ; Return
  if f1 LT f2 then return, x1 else return, x2

end
          
      
