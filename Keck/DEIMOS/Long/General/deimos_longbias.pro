;+ 
; NAME:
; deimos_longbias   
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
;   deimos_longbias, 'lblue1010.fits', ovimg
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Dec-2005 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro deimos_longbias, img, indx, ovimg, LONG=long, FITS=fits, FULL=full
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'deimos_longbias, img, ovimg, /LONG, /FULL, FITS=, [v1.1]'
      return
  endif 
  
  ;; Gain
;  [1.24, 1.20, 1.27, 1.27, 1.26, 1.25, 1.25, 1.25]

  
  ovimg = fltarr(881,4096L)
  case indx of 
      3: begin
          ;; CCD 3  ;  gain = 1.27
          fitstr = x_setfitstrct()
          ov3 = djs_median(img[2070:2130,*],1)
          f3 = x_fitrej(findgen(4096L), ov3, fitstr=fitstr)
          ovimg = img[430:1310,*] -  replicate(1.,881L)#f3
      end
      7: begin
          ;; CCD 7  ;  gain = 1.27
          ovimg = fltarr(861,4096L)
          fitstr = x_setfitstrct()
          ov7 = djs_median(img[2070:2130,*],1)
          f7 = x_fitrej(findgen(4096L), ov7, fitstr=fitstr)
          ovimg = img[770:1650,*] -  replicate(1.,881L)#f7
      end
      else: stop
  endcase

  ;; FITS
  if keyword_set( FITS ) then $
    mwrfits, ovimg, fits, /create, /silent


  return
end
