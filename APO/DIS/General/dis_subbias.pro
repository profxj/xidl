;+ 
; NAME:
; dis_subbias   
;     Version 1.1
;
; PURPOSE:
;    Subtract the overscan region from an LRIS image (blue side)
;
; CALLING SEQUENCE:
;  lrisb_subbias, img, ovimg, /LONG, FITS=, /FULL
;
; INPUTS:
;   img   -- Name of image to ov subtract
;   ovimg -- OV subtracted image
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /LONG  -- Long slit mode
;  /FULL  -- Full readout
;
; OPTIONAL OUTPUTS:
;  FITS  -- Name of fits file to write finimg to
;
; COMMENTS:
;  Currently only good for 1x1 binning
;
; EXAMPLES:
;   lrisb_subbias, 'lblue1010.fits', ovimg
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-May-2006 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro dis_subbias, img, ovimg, FITS=fits, BLUE=blue
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'dis_subbias, img, ovimg, /LONG, /FULL, FITS=, [v1.1]'
      return
  endif 
  
  if not keyword_set( BLUE ) and not keyword_set( FULL ) then begin
      ovimg = fltarr(700,2048L)
      ;; Amp3
      fitstr = x_setfitstrct()
      ov = djs_median(img[2055:*,*],1)
      fd = x_fitrej(findgen(2048L), ov, fitstr=fitstr)
      ovimg = img[0:2047L,150:849] -  replicate(1.,2048L)#fd
      ;; FITS
      if keyword_set( FITS ) then $
        mwrfits, ovimg, fits, /create, /silent
  endif


  print, 'lrisb_subbias: All done!'

  return
end
