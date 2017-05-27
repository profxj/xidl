function makearr,n,min,max,fan=fan,transpose=transpose
;+
; NAME: 
;       MAKEARR 
;
;
; PURPOSE:
;       This procedure will generate an array of lenght N which
;       runs from values MIN to MAX
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;       f = makearr(n, min, max [,fan=, transfan=, /double])
;
;
; INPUTS:
;
;       N:    The number of desired array elements in F
;       MIN:  The value of the first array element in F
;       MAX:  The value of the last array element in F
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
;       FAN:        Number of times the array is to be repeated.
;                   The final dimensions of F  will be N columns 
;                   by FAN rows.
;       /TRANSPOSE  Final dimensions of F wil be FAN columns by N 
;                   rows if FAN is specified. 
;
; OUTPUTS:
;
;       F:    Final array
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;      If you want a 5 element array which runs from 2 to 4:
;
;         IDL> f = makearr(5,2,4)
;         IDL> print, f
;             2.00000      2.50000      3.00000      3.50000      4.00000
;         
; MODIFICATION HISTORY:
; Written by John "JohnJohn" Johnson somewhere around Oct-2001
; 20 Feb 2002 JohnJohn- Added /FAN and /TRANSPOSE keywords.
; 23 Feb 2002 JohnJohn- Calculations performed in double precision. 
;                       Output in double precision if all input 
;                       parameters are double.
; 01 Mar 2002 Tim Robishaw- Spiffed up with a little Tim. 
; 08 Mar 2002 Carl suggested a change in the order of operations that
;             keeps the last number of the array equal to MAX with no
;             error.
;-

if n_params() lt 3 then begin 
    message,'Syntax: f = makearr(nelements,min,max [,fan] [,/transpose])',/info
    retall
endif

;if any of the input parameters are double, the return the answer in
;double precision.
doub = (size(n,/type) eq 5) and (size(min,/type) eq 5) and (size(max,/type) eq 5)

a = (dindgen(n)/(n-1))*(max-min)+min

if n_elements(fan) ne 0 then a = a##(dblarr(fan)+1)
if KEYWORD_SET(transpose) then a = transpose(a)

if not doub then a = float(a)
return,a
end
