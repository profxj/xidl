;+
; NAME:
;   long_slits2mask
;
; PURPOSE:
;   Convert trace sets that describe slit positions into a mask
;   with the slit number for each pixel
;
; CALLING SEQUENCE:
;   slitmask = long_slits2mask(tset_slits, [ xshift=, nslit=
;                                        , onlyslits =] )
;
; INPUTS:
;   tset_slits - 2-element array of trace sets, where the first defines
;                the starting slit positions, and the second one defines
;                the ending slit positions
;
; OPTIONAL INPUTS:
;   xshift     - Number of pixels to shift the traces in the output image;
;                this can bbe a fractional number
;   onlyslits  - only return a mask with specified slit numbers
;
; OUTPUTS:
;   slitmask   - Mask image, with values of zero where there is no slit,
;                and the slit number (starting at 1) for each object
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
;
; REVISION HISTORY:
;   10-Mar-2005  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
function long_slits2mask, tset_slits1, xshift=xshift, nslit=nslit $
                          , VERBOSE = VERBOSE, ONLYSLITS = ONLYSLITS
  
   if (size(tset_slits1,/tname) NE 'STRUCT' $
    OR n_elements(tset_slits1) NE 2) then $
    message, 'TSET_SLITS must be a 2-element structure'

   tset_slits = tset_slits1
   if (keyword_set(xshift)) then $    
    tset_slits.coeff[0,*] = tset_slits.coeff[0,*] + xshift

   dims = tset_slits[0].dims
   nx = dims[0]
   ny = dims[1]

   traceset2xy, tset_slits[0], yy1, xx1
   traceset2xy, tset_slits[1], yy2, xx2
   
   if (size(xx1,/n_dimen) EQ 1) then nslit = 1 $
    else nslit = (size(xx1,/dimens))[1]

   IF KEYWORD_SET(ONLYSLITS) THEN BEGIN
      nloop = n_elements(onlyslits)
      doslits = onlyslits - 1L
   ENDIF ELSE BEGIN
      nloop = nslit
      doslits = lindgen(nslit)
   ENDELSE

   ; Generate the mask image
   slitmask = intarr(dims)
   for ii=0L, nloop-1L do begin
      islit = doslits[ii] 
      for iy=0L, ny-1L do begin
         x1 = round(xx1[iy,islit])
         x2 = round(xx2[iy,islit])
         if (x1 GE x2) AND KEYWORD_SET(VERBOSE) THEN $
          splog, 'WARNING: Slit start and end positions appear to cross!'
         x1 = x1 > 0
         x2 = x2 < (nx-1)
         if (x1 LE x2) then begin
             if (total(slitmask[x1:x2, iy]) GT 0) $
               AND KEYWORD_SET(VERBOSE) then $
               splog, 'WARNING: Slits appear to overlap!'
             slitmask[x1:x2, iy] = islit+1
         endif
      endfor
   endfor

   return, slitmask
end
;------------------------------------------------------------------------------
