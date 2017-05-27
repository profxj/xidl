;+
; NAME:
;   deimos_rectify_slit
;
; PURPOSE:
;   extract a rectified 2d spectrum from 1-chip DEIMOS image,
;   transpose result so that first axis is spectral direction
;
; CALLING SEQUENCE:
;   spec = deimos_rectify_slit(im, ivar, x0, x1, recenter=recenter, $
;             npad=npad, interp=interp, xshift=xshift, mask=mask)
;
; INPUTS:
;   im         - 2D spectral image
;   ivar       - inverse variance of image
;   x0         - "left" side of region (array [ny])
;   x1         - "right" side of region (array [ny])
;
; OPTIONAL KEYWORDS:
;   npad       - pad extracted region with npad pixels on each side
;   recenter   - attempt to recenter the spectrum (flux-weighted)
;   interp     - sample fractional pixel locations with
;                lin. interpolation
; OUTPUTS:
;   spec       - array[nx+2*npad, ny] with extracted region
;
; OPTIONAL OUTPUTS:
;   xshift     - amount shifted in the recenter operation
;   mask       - same size as spec.  1=good, 0=pixel out of image.
;
; COMMENTS:
;   An inverse variance array should be passed, but if it is not the
;   procedure pretends it was all 1.  
;
; BUGS:
;   /recenter does not work very well if the trace falls off the edge
;   of the image substantially. 
;
; EXAMPLES:
;
; REVISION HISTORY:
;   01-Dec-2000  Written by D. Finkbeiner, Berkeley
;   11-Apr-2002  transpose output, MD
;-
;------------------------------------------------------------------------------
function deimos_rectify_slit, im, ivar, x0, x1, npad=npad, $
             recenter=recenter, interp=interp, xshift=xshift, mask=mask

  sizex = (size(im, /dimens))[0]
  if not keyword_set(npad) then npad = 0 ;extract npad extra pix on each side

  ny = n_elements(x0)
  if (size(im, /dimens))[1] NE ny then message, 'wrong size!'

  nx = ceil(median(abs(x1-x0)))+npad*2  ; median width
  if nx LE npad*2 then message, 'invalid!'

  kernel=9+2*(nx/300)

  xcorrs=findgen(2*kernel-5)-kernel+3
  xcorrs=xcorrs*2.5
  if nx lt 200 then xcorrs=xcorrs[2*findgen(kernel-2)]

; recenter
  if keyword_set(recenter) then begin 
     ; call recursively
     specpad = deimos_rectify_slit(im, ivar, x0, x1, $
                                   npad=(npad+2) > 4)
     wid = (size(specpad, /dim))[1]
     mid = specpad[*, wid/2]
     smid = smooth(mid, kernel)
;     smid=dilate(mid gt median(mid)+stdev(mid),fltarr(2*kernel+1)+1)


     cc = fltarr(wid)

     for i=0, wid-1 do cc[i] = max(c_correlate(smid, specpad[*, i], xcorrs))

     cc2=total(specpad,1)
     cc2m=djs_median(cc2,width=41<(wid/10),boundary='reflect')>1.E-5
     
     dx = findgen(wid)+(0.5-wid/2.)

; DF's code
     xshift = total(dx*cc)/total(cc-min(cc))

; but that's not a proper centroid, so try:
     xshift2 = total(dx*(cc-min(cc)))/total(cc-min(cc))

; the following _might_ work better on long slits
     xshift3 = total(dx*cc2/cc2m)/total(cc2/cc2m)

;     print,xshift,xshift2,xshift3

     xshift=xshift2

  endif else xshift = 0

  xav  = (x0+x1)/2.
  spec = fltarr(nx, ny)
  ind0 = xav+(0.5+xshift-nx/2.)
  indx = lindgen(nx)

; if a mask is desired:
  if arg_present(mask) then begin
     mask = bytarr(nx, ny)
     if keyword_set(interp) then begin 
        for i=0, ny-1 do mask[*, i] = interpolate(float(ivar[*, i] NE 0), ind0[i]+indx) EQ 1
     endif else begin 
        for i=0, ny-1 do mask[*, i] = ivar[round(ind0[i])+indx, i] NE 0
     endelse 
     mask = transpose(mask)
  endif

; extract region
  if keyword_set(interp) then begin 
     for i=0, ny-1 do spec[*, i] = interpolate(im[*, i], ind0[i]+indx)
  endif else begin 
     for i=0, ny-1 do spec[*, i] = im[round(ind0[i])+indx, i]
  endelse 


  return, transpose(spec)
end
