;+ 
; NAME:
; x_imarith   
;    Version 1.1
;
; PURPOSE:
;    Performs arithmetic on two fits images (akin to IRAF)
;
; CALLING SEQUENCE:
;   x_imarith, img1, oper, img2, outimg, FITS=
;
; INPUTS:
;   img1 -- Image 1
;   img2 -- Image 2
;   oper -- String math operation (+,*,/,-)
;
; RETURNS:
;
; OUTPUTS:
;   fimg -- Output image 
;
; OPTIONAL KEYWORDS:
;
;
; OPTIONAL OUTPUTS:
;   /FITS - fimg is a fits file (string)
;
; COMMENTS:
;
; EXAMPLES:
;   x_imarith, 'f1.fits', '+', 'f2.fits', fimg
;   x_imarith, 'f1.fits', '+', 'f2.fits', 's12.fits', /FITS
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   03-Aug-2001 Written by JXP
;   29-Dec-2001 Modified
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_imarith, img1, oper, img2, fimg, FITS=fits

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'x_imarith, img1, oper, img2, fimg, /FITS (v1.1)'
      return
  endif 
  
  i1 = x_readimg(img1, HEAD=h1, /fscale)
  i2 = x_readimg(img2, /fscale)


  case oper of 
      '+' : outimg = i1 + i2
      '-' : outimg = i1 - i2
      '/' : begin
          outimg = i1 * 0.
          gd = where(i2 NE 0., ngd)
          if ngd NE 0 then outimg[gd] = i1[gd] / i2[gd]
      end 
      '*' : outimg = i1 * i2
      else : begin
          print, 'Operation', oper, ' not defined'
          return
      end
  endcase


; FITS file

  if keyword_set( FITS ) then begin
      sxaddhist, strjoin(['x_imarith: ', img1, oper, img2]), h1
      mwrfits, outimg, fimg, h1, /silent, /create
  endif else fimg = temporary(outimg)

end
