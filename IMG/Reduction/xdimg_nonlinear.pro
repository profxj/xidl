;+ 
; NAME:
; xdimg_nonlinear   
;  Version 1.0
;
; PURPOSE:
;    Performs non-linear correction on a ccd image.  Only setup
;   for SITe3 at this point.
;
; CALLING SEQUENCE:
;   newimg = xdimg_nonlinear(img, ccd)
;
; INPUTS:
;    img  -  image
;    ccd  -  (SITe3)
;
; RETURNS:
;    newimg - Corrected image
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
;   newimg = xdimg_nonlinear(img, 'SITe3')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function xdimg_nonlinear, img, ccd

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'newimg = xdimg_nonlinear(img, ccd) (v1.1)'
      return, -1
  endif 
  
;  Optional Keywords

;  CCD dependent

  sz = size(img, /dimensions)
  case ccd of
      'SITe3': begin
          newimg = fltarr(sz[0], sz[1])
          a = where(img GT 1.)
          newimg[a] = (1 + (0.127*alog10(img[a]))^4)*img[a]
      end
      else: return, img
  endcase

  return, newimg
end

