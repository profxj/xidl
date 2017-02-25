;+ 
; NAME:
; x_ordrectify
;     Version 1.1
;
; PURPOSE:
;  Return an image of a single order, which is rectified,
;  has a width of 2*long(min(rhedg-lhedg)/2) + 1
;  and conserves counts in each row.
;  lhedg maps to the first column, rhedg to last column
;  and center of order to central column
;
; CALLING SEQUENCE:
;   
;  rect_image = x_ordrectify(arc_img, lhedg, rhedg, /NOCORRECT) 
;
; INPUTS:
;   img  -  Raw image
;   lhedg - Left hand edge of the order
;   rhedg - Right hand edge of the order
;
; RETURNS:
;   rect_img  -  Rectified 2D image
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOCORRECT  - ??
;  HALFW - Half width of the order (default: Minimum half width)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  rect_image = x_ordrectify(arc_img, lhedg, rhedg, /NOCORRECT)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_ordrectify
;
; REVISION HISTORY:
;   ??--2004 Written by SB
;-
;------------------------------------------------------------------------------
function x_ordrectify, img, lhedg, rhedg, nocorrect=nocorrect, $
       halfw=halfw

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rect_image = x_ordrectify(arc_img, lhedg, rhedg, /NOCORRECT) [v1.1] ' 
      return, -1
  endif

;  nord = n_elements(ordr_str)

  sz_img = size(img, /dimen)
  
  width = rhedg - lhedg
  if NOT keyword_set(halfw) then halfw = min(width)/2.

  girth = 2L*long(halfw) + 1L
  x =  findgen(girth)/(girth-1)

  slit_pos = x # width + replicate(1.0,girth) # (lhedg)
  left_indx = long(slit_pos)
  right_indx = left_indx + 1
  row_indx   = replicate(1,girth) # lindgen(sz_img[1])
  right_frac = (slit_pos - left_indx)
  left_frac = 1.0 - right_frac
  
  ;; Correction
  if keyword_set(nocorrect) then correction=1.0 $
  else correction = replicate(1./girth,girth) # width
  left_frac = left_frac * correction 
  right_frac = right_frac * correction 
  
  rect_image = fltarr(girth, sz_img[1])
  good = where(left_indx GT 0 AND right_indx LT sz_img[0])

  ;; Interpolate over all pixels in the rectified image
  if good[0] NE -1 then begin
      rect_image[good] = left_frac[good]*img[left_indx[good],row_indx[good]]+$
        right_frac[good]*img[right_indx[good],row_indx[good]]
  endif
  
  return, rect_image
end                   
    
