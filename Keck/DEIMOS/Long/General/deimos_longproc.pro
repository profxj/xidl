;+ 
; NAME:
; deimos_longproc
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
pro deimos_longproc, img, indx, finimg, FITS=fits, FLAT=flat


  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'deimos_longproc, img, indx, [finimg], FLAT=, FITS=  [v1.1]'
      return
  endif 

  ;; Gain
  gain = [1.24, 1.20, 1.27, 1.27, 1.26, 1.25, 1.25, 1.25]
  
  ;; IMG 3 FIRST ;;

  ;; Open image
  print, 'deimos_longproc: Reading...'
  raw = xmrdfits(img, indx, /fscale, /silent)
  
; Bias subtract
  print, 'deimos_longproc: Bias subtracting'
  deimos_longbias, raw, indx, ovimg

; Add in the gain
  ovimg = ovimg * gain[indx-1]

; Flatten
  print, 'deimos_longproc: Flattening'
  if keyword_set( FLAT ) then begin
      fdat = xmrdfits(flat, /silent)
      finimg = ovimg / fdat
  endif else finimg = ovimg

; Output
  if keyword_set( FITS ) then mwrfits, finimg, fits, /create

  print, 'deimos_longproc: All done!'

  return
end
