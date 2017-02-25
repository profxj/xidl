;+
; NAME:
;   x_medsigclip
;
; PURPOSE:
;   Median multiple images with sigma-rejection.  Akin to avsigclip
;   but does median analysis instead of averaging.  This is my favorite
;   combine routine.
;
; CALLING SEQUENCE:
;   result = x_medsigclip( array, [ dimension, siglo=, sighi=, maxiter=, $
;    inmask=, ] )
;
; INPUTS:
;   array      - N-dimensional array
;
; OPTIONAL INPUTS:
;   dimension  - The dimension over which to collapse the data.  If ommitted,
;                assume that the last dimension is the one to collapse.
;   siglo     - Low Sigma for rejection; default to 3.0.
;   sighi     - High Sigma for rejection; default to 3.0.
;   maxiter    - Maximum number of sigma rejection iterations.  One iteration
;                does no sigma rejection; default to 10 iterations.
;   inmask     - Input mask, setting =0 for good elements
;   gain       - Gain; Crucial for calculating sigma; default to 1.0
;   rn         -  Readnoise; Crucial for calculating sigma; default to 7.0
;
; OUTPUTS:
;   result     - The output array.
;   outmask    - Output mask, setting =0 for good elements, =1 for bad.
;                Any pixels masked in INMASK are also masked in OUTMASK.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The DIMENSION input is analogous to that used by the IDL built-in
;   function TOTAL.
;
; EXAMPLES:
;   Create a data cube of 10 random-valued 100x200 images.  At each pixel in
;   the image, compute the median of the 10 values, but rejecting 3-sigma
;   outliers:
;   > array = randomu(123,100,200,10)
;   > med = x_medsigclip(array, siglo=3., sighi=4.)
;
;
;   If all points are masked in any given vector or array, a mean and
;   dispersion are computed for all the points.  Is this the behaviour we want?
;   If you want to replace those values with zeros instead, look at OUTMASK:
;   > array = randomu(123,100,200)
;   > inmask = bytarr(100,200)
;   > inmask[*,8] = 1 ; mask all of row #8
;   > ave = x_avsigclip(array, 1, inmask=inmask, outmask=outmask)
;   > ibad = where( total(1-outmask, 1) EQ 0)
;   > if (ibad[0] NE -1) then ave[ibad] = 0 ; zero-out bad rows
;
; BUGS:
;
; PROCEDURES CALLED:
;   Dynamic link to arravsigclip.c
;
; REVISION HISTORY:
;   07-Jul-1999  Based on djs_avsigclip; Written by David Schlegel, Princeton.
;   28-Jul-2001  Modfied by JXP (Added lo/hi sigrej, does median)
;-
;------------------------------------------------------------------------------
function x_medsigclip, array, dim, siglo=siglo, sighi=sighi, maxiter=maxiter, $
 INMASK=inmask, GAIN=gain, RN=rn

   ; Need at least 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - result = x_medsigclip( array, [ dimension, sigrej=, maxiter=, $'
      print, ' inmask= ])';
      return, -1
   endif

   if (NOT keyword_set(dim)) then dim = size(array, /n_dim)
   if (NOT keyword_set(sighi)) then sighi = 3.0
   if (NOT keyword_set(siglo)) then siglo = 3.0
   if (NOT keyword_set(gain)) then gain = 1.0
   if (NOT keyword_set(rn)) then rn = 7.0
   if (NOT keyword_set(maxiter)) then maxiter = 5

   sz = N_elements(array)
   dimvec = size(array, /dimensions)
   ndim = N_elements(dimvec)

   if (dim GT ndim OR dim LT 1) then begin
      message, 'DIM must be between 1 and '+string(ndim)+' inclusive'
   endif

   ; Allocate memory for the output array
   if (ndim GT 1) then $
    newdimvec = dimvec[ where(lindgen(ndim)+1 NE dim) ] $
   else $
    newdimvec = [1]
   newsize = N_elements(array) / dimvec[dim-1]
   medarr = reform(fltarr(newsize), newdimvec)

   soname = filepath('libxmath.' + idlutils_so_ext(), $
    root_dir=getenv('XIDL_DIR'), subdirectory='/lib')
   soname2 = filepath('libxmath.' + idlutils_so_ext(), $
    root_dir=getenv('XIDL_DIR'), subdirectory='/lib')

   ; Check data types.  Will be strict in requiring identical types!
   if keyword_set(inmask) then begin

       if( size(array, /type) NE 4  OR $
           size(inmask, /type) NE 1 ) then begin
           print, 'x_medsigclip: Data types not right!'
       endif
       retval = call_external(soname2, 'arrmedsigmask', $
                              ndim, dimvec, array, long(dim), float(siglo), $
                              float(sighi), long(maxiter), float(gain), float(rn), $
                              medarr, inmask)

   endif else begin

       if size(array, /type) NE 4  then begin
           print, 'x_medsigclip: Data types not right!'
       endif
       retval = call_external(soname, 'arrmedsigclip', $
                              ndim, dimvec, array, long(dim), float(siglo), $
                              float(sighi), long(maxiter), float(gain), float(rn), $
                              medarr)
   endelse

   return, medarr
end
;------------------------------------------------------------------------------
