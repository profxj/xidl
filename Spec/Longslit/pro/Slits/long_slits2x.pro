;+
; NAME:
;   long_slits2x
;
; PURPOSE:
;   Find the location of objects within each slit mask
;
; CALLING SEQUENCE:
;   ximg = long_slits2x( tset_slits, [ slitid=, xshift=, nslit= ] )
;
; INPUTS:
;   tset_slits - Trace sets with slit start/end positions
;   slitid     - Slit ID number(s) (1-indexed); default to all slit IDs
;   xshift     - Number of pixels to shift the traces in the output image;
;                this can bbe a fractional number
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   ximg       - Image with the fractional X (spatial) position along
;                the specified slit for each pixel
;
; OPTIONAL OUTPUTS:
;   nslit      - Number of slits
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   traceset2xy
;   
; REVISION HISTORY:
;   11-Mar-2005  Written by D. Schlegel, LBL
;-  
;------------------------------------------------------------------------------
function long_slits2x, tset_slits1, pixleft = pixleft, pixright = pixright $
                       , slitid = slitid1, TOL_EDG = TOL_EDG1 $
                       , EDGMASK = EDGMASK $
                       , xshift = xshift, nslit = nslit

   if (size(tset_slits1,/tname) NE 'STRUCT' $
    OR n_elements(tset_slits1) NE 2) then $
    message, 'TSET_SLITS must be a 2-element structure'
   IF n_elements(TOL_EDG1) EQ 0 THEN TOL_EDG = [3L, 3L] $
   ELSE BEGIN
      IF n_elements(tol_edg1) EQ 1 THEN tol_edg = [tol_edg1, tol_edg1] $
      ELSE IF n_elements(tol_edg1) EQ 2 THEN tol_edg = tol_edg1 $
      ELSE message, 'incorrect number of elementns for tol_edg'
   ENDELSE
   
   tset_slits = tset_slits1
   if (keyword_set(xshift)) then $
    tset_slits.coeff[0,*] = tset_slits.coeff[0,*] + xshift

   traceset2xy, tset_slits[0], yy1, xx1
   traceset2xy, tset_slits[1], yy2, xx2

   ;; Generate the output image
   dims = tset_slits[0].dims
   ximg = fltarr(dims)
   pixleft = fltarr(dims)
   pixright = fltarr(dims)
   if (size(xx1,/n_dimen) EQ 1) then nslit = 1 $
    else nslit = (size(xx1,/dimens))[1]

   if (keyword_set(slitid1)) then doslits = slitid1 $
    else doslits = lindgen(nslit) + 1

   ;; Loop over each slit
   for islit = 0L, n_elements(doslits)-1 do begin
       slitid = doslits[islit]
       ;; How many pixels wide is the slit at each Y?
       xsize = xx2[*, slitid-1] - xx1[*, slitid-1]
       badp = where(xsize LE 0., nbad)
       if nbad GT 0 then begin
           meds = median(xsize)
           print, 'long_slits2x:  Something goofy in slit #'+strtrim(slitid,2)
           print, 'long_slits2x:  Probably a bad slit (e.g. a star box)'
           print, 'long_slits2x:  It is best to expunge this slit with long_fiddleslits'
           print, 'long_slits2x:  Proceed at your own risk, with a slit width of ', meds
           print, 'long_slits2x:  Or set meds to your liking'
           stop
           xx2[*,slitid-1] = xx1[*,slitid-1] + meds
       endif
       
       for iy = 0L, dims[1]-1L do begin
           ix1 = ceil(xx1[iy, slitid-1]) > 0
           ix2 = floor(xx2[iy, slitid-1]) < (dims[0]-1)
           ix1 = ix1 < (dims[0]-1)
           ix2 = ix2 > 0
           ximg[ix1:ix2, iy] = (findgen(ix2-ix1+1) + ix1 - xx1[iy, slitid-1]) $
             / xsize[iy]
           pixleft[ix1:ix2, iy] = (findgen(ix2-ix1+1) + ix1 - xx1[iy, slitid-1])
           pixright[ix1:ix2, iy] = (xx2[iy, slitid-1] - ix2 + $
                                    reverse(findgen(ix2-ix1+1)))
       endfor
   endfor
   slitmask = long_slits2mask(tset_slits)
   edgmask  = slitmask GT 0 AND (pixleft LT tol_edg[0] OR $
                                 pixright LT tol_edg[1])
   return, ximg
end
;------------------------------------------------------------------------------
