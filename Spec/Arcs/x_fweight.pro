;+
; NAME:
;   x_fweight
;
; PURPOSE:
;   Recenter a x using flux-weighting centers.
;
; CALLING SEQUENCE:
;   xnew = x_fweight( fimage, xcen, ycen, [radius=radius, xerr=xerr, 
;    invvar=invvar] )
;
; INPUTS:
;   fimage     - Image
;   xcen       - Initial guesses for X centers
;   ycen       - Y positions corresponding to "xcen" (expected as integers)
;
; OPTIONAL KEYWORDS:
;   radius     - Radius for centroiding; default to 3.0
;   invvar     - Inverse variance of image used only in computing errors XERR.
;                If not set, then INVVAR=1 is used.
;
; OUTPUTS:
;   xnew       - New X centers
;
; OPTIONAL OUTPUTS:
;   xerr       - Formal errors for XNEW; set equal to 999.0 if there are any
;                masked pixels in a centroiding region (e.g., if INVVAR=0)
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   Dynamic link to x_fweight.c
;
; REVISION HISTORY:
;   24-Mar-1999  Written by David Schlegel, Princeton.
;   14-Aug-1999  Revised to loop
;-
;------------------------------------------------------------------------------
function x_fweight, fimage, xcen, ycen, radius=radius, xerr=xerr, $
 invvar=invvar, NLOOP=nloop

   ; Need 3 parameters
   if (N_params() LT 3) then begin
      print, 'Syntax - xnew = x_fweight( fimage, xcen, ycen, [radius=, $'
      print, ' xerr=, invvar=, NLOOP=] )'
      return, -1
   endif
   if (NOT keyword_set(radius)) then radius = 3.0
   if NOT keyword_set(radius) then nloop = 1L

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   ncen = N_elements(xcen)
   xnew = float(xcen)
   xerr = 0.0 * xnew ; Allocate memory

   if (NOT keyword_set(invvar)) then invvar = 0.0 * fimage + 1.0

   soname = filepath('libspec2d.so', $
                     root_dir=getenv('IDLSPEC2D_DIR'), subdirectory='lib')
   result = call_external(soname, 'trace_xfweight', $
                          nx, ny, float(fimage), float(invvar), $
                          float(radius), ncen, xnew, long(ycen), xerr, $
                          nloop)
       
   return, xnew
end
;------------------------------------------------------------------------------
