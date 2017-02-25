;+ 
; NAME:
; dis_proc
;     Version 1.1
;
; PURPOSE:
;    Process the image (bias subtract, flatten and multiply by gain)
;
; CALLING SEQUENCE:
;  lrisb_proc, img, flat, finimg, /LONG, FITS=, /FULL
;
; INPUTS:
;   img  -- Name of file to process
;   flat -- Name of flat
;
; RETURNS:
;  finimg -- Processed image file
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /LONG  -- Long slit mode
;  /FULL  -- Keyword for lrisb_subbias
;
; OPTIONAL OUTPUTS:
;  FITS  -- Name of fits file to write finimg to
;
; COMMENTS:
;  Currently only good for 1x1 binning
;
; EXAMPLES:
;   lrisb_subbias, lris, [47L,48L,49L]
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-May-2006 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro dis_proc, img, finimg, flat, FITS=fits,  NOFLAT=noflat, BLUE=blue


  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'dis_proc, img, finimg, [flat], /LONG, FITS=, /BLUE [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( GAIN ) then begin
      if keyword_set(BLUE) then gain = 1.74 else gain = 1.89
  endif

; Open image
  print, 'dis_proc: Reading...'
  raw = xmrdfits(img, /fscale, /silent)
  
; Bias subtract
  print, 'dis_proc: Bias subtracting'
  dis_subbias, raw, ovimg, BLUE=blue

; Add in the gain
  ovimg = ovimg * gain

; Flatten
  print, 'dis_proc: Flattening'
  if not keyword_set( NOFLAT ) then begin
      fdat = xmrdfits(flat, /silent)
      finimg = ovimg / fdat
  endif else finimg = ovimg

; Output
  if keyword_set( FITS ) then mwrfits, finimg, fits, /create

  print, 'dis_proc: All done!'

  return
end
