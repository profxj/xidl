;+
; NAME:
;       skysub
;
; PURPOSE:
;       subtract the sky from a 2D spectrum using the b_spline method
;
; CALLING SEQUENCE:
;       skysub, obj, skyarr, Lxyo, Lxys, skyrows=skyrows, objnosky=objnosky
;
; INPUTS:
;       obj     - 2D normalized spectrum (returned by procedure norma) 
;       skyarr  - 2D array containing the sky 
;       Lxyo    - wavelength solution for the 2D spectrum Lambda=lambda(xy)
;       Lxys    - wavelength solution for the sky array
;	invvar 	- inverse variance of the data, (multiplied by mask) 
;
;
; OPTIONAL INPUTS:
;      
;
; REQUIRED KEYWORDS:
;      objnosky  - sky subtracted 2D spectrum
;
; OPTIONAL KEYWORDS: 
;      skyrows   - list of subscripts specifying the rows of skyarr to
;                  be used for creating the bsplined sky 
;   
;
; OUTPUTS:
;   
;
; OPTIONAL OUTPUTS:
;   
;   
; COMMENTS:
;
;
; EXAMPLES:
;
;
; BUGS:
;
;
; PROCEDURES CALLED:
;   
;
; REVISION HISTORY:
;    19-Mar-2001 Written by Chris Marinoni, Berkeley
;    30-Aug_20001 modified by Andrew Sheinis, UCSC to accept ESI data, 4096 columns
;    30-apr_20001 modified by Andrew Sheinis, added rejection.	   
;-
;--------------------------------------------------------------------------- 

pro skysub_esi, obj, skyarr, Lxyo, Lxys, skyrows=skyrows, objnosky=objnosky, range=range,invvar=invvar


   if Not keyword_set(skyrows) then skyrows = findgen((size(skyarr))(2))
  
   if max(skyrows) gt (size(skyarr))(2)-1 or min(skyrows) lt 0  then begin
     print, 'ERROR: skyrow subscript out of range'
     stop
   endif
  
   s = skyrows
   lxysext = lxys(*, s)             
   isext   = skyarr(*, s)

   ; data ordering  
   isort = bsort(lxysext)
   l1d = lxysext(isort)
   i1d = isext(isort)

   ; define a vector of breakpoints 
   bkpt =  lxys(1:4095, 20)      ;THIS IS NOT RIGHT! ais 
   everyn=.667*n_elements(skyrows)

   a = systime(1)
;   error_code = bspline_fit(l1d, i1d, 1./sqrt(i1d), bset, fullbkpt = bkpt)

bset = bspline_iterfit(l1d, i1d,everyn=everyn,maxiter=3,invvar=invvar,upper=3,lower=5,yfit=yfit)
;bset = bspline_iterfit(l1d, i1d,everyn=everyn,maxiter=3,invvar=invvar,upper=3,lower=3,yfit=yfit)

   skyobj = bspline_valu(lxyo, bset)  

    ; sky subtraction
   objnosky = obj-skyobj
   print,  'calculation time', systime(1)-a

;stop

return
end
      




