;+ 
; NAME:
; x_addtwoflats
;   Version 1.1
;
; PURPOSE:
;    Combines two flats, rejecting Cosmic Rays.  This routine is
;   used extensively in the WFCCD routines.
;
; CALLING SEQUENCE:
;   
;   flat = x_addtwoflats, img1, img2, VAR=, GAIN=, RN=
;
; INPUTS:
;   img1       - OV subtracted img1 (data or fits)
;   img2       - OV subtracted img2 (data or fits)
;
; RETURNS:
;
; OUTPUTS:
;   flat      - Combined flat with CR rejected
;
; OPTIONAL KEYWORDS:
;  GAIN       - gain
;  RN         - Read noise (for the VAR output only)
;
; OPTIONAL OUTPUTS:
;   VAR       - Variance image
;
; COMMENTS:
;
; EXAMPLES:
;   flat = x_addtwoflats(img1, img2)
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Jan-2002 Written by JXP
;-
;------------------------------------------------------------------------------

function x_addtwoflats, img1, img2, GAIN=gain, RN=rn, VAR=var

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'flat = x_addtwoflats(img1, img2, GAIN=) [v1.0]'
    return, -1
  endif 


;  Optional Keywords

  if not keyword_set( RN )   then rn = 5.
  if not keyword_set( GAIN ) then gain = 1.
  if not keyword_set( SATUR )then satur = 30000.
  

; Allow img to be fits file or data

  dat1 = x_readimg(img1, /fscale)
  dat2 = x_readimg(img2, /fscale)

  ; Ratio
  rtio = dat1/dat2
  ; Stats on the ratio
  djs_iterstat, rtio, sigrej=3.0, median=med_rtio, sigma=sig_rtio, maxiter=2

  ; Find all bad pixels
  bdpix = where(abs(rtio-med_rtio) GT 5.*sig_rtio, nbad)

  ; Average (not add)
  flat = (dat1+dat2)/2.

  ; Variance
  if arg_present(VAR) then var = ( flat/2. + rn^2)

  ; Take minimum of bad pixels
  if nbad GT 0 then begin
      flat[bdpix] = dat1[bdpix] < dat2[bdpix]
      if arg_present(VAR) then var[bdpix] = ( flat[bdpix] + rn^2 )
  endif

  delvarx, dat1, dat2, rtio, bdpix


  return, flat

end
