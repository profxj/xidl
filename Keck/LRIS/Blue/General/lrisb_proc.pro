;+ 
; NAME:
; lrisb_proc
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
;   30-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lrisb_proc, img, flat, finimg, LONG=long, FITS=fits, FULL=full, $
                NOFLAT=noflat


  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'lrisb_proc, img, flat, finimg, /LONG, FITS=, /FULL [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( GAIN ) then gain = 1.54
;  if not keyword_set( OVROOT ) then ovroot = 'OV/ov_'
;  if not keyword_set( BIASFIL ) then begin
;      if lrisb[indx[0]].mode EQ 0 then biasfil = 'Bias/BiasI.fits' $
;        else biasfil = 'Bias/BiasS.fits'
;  endif

; Open image
  if not keyword_set( LONG ) then stop
  print, 'lrisb_proc: Reading...'
  raw = xmrdfits(img, /fscale, /silent)
  
; Bias subtract
  print, 'lrisb_proc: Bias subtracting'
  lrisb_subbias, raw, ovimg, LONG=long, FULL=full

; Add in the gain
  ovimg = ovimg * gain

; Flatten
  print, 'lrisb_proc: Flattening'
  if not keyword_set( NOFLAT ) then begin
      fdat = xmrdfits(flat, /silent)
      finimg = ovimg / fdat
  endif else finimg = ovimg

; Output
  if keyword_set( FITS ) then mwrfits, finimg, fits, /create

  print, 'lrisb_proc: All done!'

  return
end
