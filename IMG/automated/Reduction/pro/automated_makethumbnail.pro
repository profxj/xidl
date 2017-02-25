;+
; NAME:
;   automated_makethumbnail
;
; PURPOSE:
; 
; Program outputs a simple thumbnail jpg image for later inspection or 
; for putting up on a website
; 
;
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; 
; OPTIONAL INPUTS:
;  
; 
; OUTPUTS: 
;
;
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
; 
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;           2011  Written by KK
;  
;-
;------------------------------------------------------------------------------
;; 



PRO make_thumbnail, image, output
  
  z_scale = ZSCALE_RANGE(image)                 ;Grab z-scale min and max
  image = BYTSCL(image, z_scale[0], z_scale[1]) ;Scale image pixel values so it's easy to see
  image = REBIN(image, 512, 512)                ;Resize thumbnail image to 512x512
  image = REVERSE(image, 1)                     ;Flip image horizontally so that it matches 
  
  WRITE_JPEG, output, image
  
END
