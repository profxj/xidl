;+ 
; NAME:
; lrisb_subbias   
;     Version 1.1
;
; PURPOSE:
;    Median combine all ZRO frames (darks)
;      WARNING!  Assumes images are all of 1 mode (e.g. IMG, ECH, LWD)!!
;
; CALLING SEQUENCE:
;   
;  lrisb_subbias, lris, indx
;
; INPUTS:
;   lrisb   -  ESI structure
;   indx  -  Index numbers of frame to subtract (default output is OV)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  BIASFIL= - Name of bias file (default: Bias/BiasS[I].fits)
;  OVROOT=  - Root name of OV file (default: OV/ov_ )
;  /FORCE   - Overwrite existing OV files 
;
; OPTIONAL OUTPUTS:
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

pro lrisb_proc, img, flat, finimg, LONG=long, FITS=fits


  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'lrisb_subbias, img, flat, finimg '
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
  lrisb_subbias, raw, ovimg, LONG=long

; Add in the gain
  ovimg = ovimg * gain

; Flatten
  print, 'lrisb_proc: Flattening'
  fdat = xmrdfits(flat, /silent)

  finimg = ovimg / fdat

; Output
  if keyword_set( FITS ) then mwrfits, finimg, fits, /create

  print, 'lrisb_proc: All done!'

  return
end
