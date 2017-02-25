;+ 
; NAME:
; x_normspec   
;    Version 1.0
;
; PURPOSE:
;    Normalizes a multi-slit image given a normalized flat
;
; CALLING SEQUENCE:
;   nrmspec = x_normspec( img, flat, [var, nrmvar])
;
; INPUTS:
;   flat       - Flat image or fits file
;   slistr     - Slit structure
;
; RETURNS:
;   nrmspec    - Normalized flat   
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   nrmspec = x_normspec( img, flat)
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_normspec, img, flat, var, nrmvar, PIX=pix


;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'nrmspec = x_normspec(img, flat, [var,nrmvar], PIX=) [v1.0]'
    return, -1
  endif 


;  Optional Keywords

; Allow flat,img to be fits file

;  dimg = x_readimg(img, /fscale)
;  dflat = x_readimg(flat, /fscale)
  sz_img = size(img, /dimensions)

; Create the output image
  nrmimg = make_array(sz_img[0], sz_img[1], type=size(img,/type))
  if keyword_set(VAR) then begin
      nrmvar = make_array(sz_img[0], sz_img[1], type=size(var,/type))
      nrmvar[*] = -1.
  endif
  
; Avoid dividing by 0
  if keyword_set( PIX ) then a=pix $
  else a = where(flat NE 0)

  nrmimg[a] = img[a]/flat[a]

; VAR
  if keyword_set( VAR ) then $
    nrmvar[a] = var[a] / (flat[a])^2

  return, nrmimg

end
  
  
