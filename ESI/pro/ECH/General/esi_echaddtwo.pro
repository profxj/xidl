;+ 
; NAME:
; esi_echaddtwo
;   Version 1.0
;
; PURPOSE:
;    Combines two flats, rejecting Cosmic Rays
;    Assumes images have nearly the same exposure time
;
; CALLING SEQUENCE:
;   
;   img = esi_echaddtwo(esi, indx, VAR=var)
;
; INPUTS:
;   esi
;   indx
;
; RETURNS:
;   img       - Combine image
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   VAR       - Variance
;
; COMMENTS:
;
; EXAMPLES:
;   img = esi_echaddtwo(esi, indx)
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   19-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------

function esi_echaddtwo, esi, indx, VAR=var, SCALE=scale

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'img = esi_echaddtwo(esi, indx, VAR=, /SCALE) [v1.0]'
    return, -1
  endif 


;  Optional Keywords

; Allow img to be fits file or data

  dat1 = xmrdfits(esi[indx[0]].img_final, 2, /silent)
  dat2 = xmrdfits(esi[indx[1]].img_final, 2, /silent)

; Variance
  if arg_present(VAR) then begin
      var1 = xmrdfits(esi[indx[0]].img_final, 1, /silent)
      var2 = xmrdfits(esi[indx[1]].img_final, 1, /silent)
  endif

  if keyword_set( SCALE ) then begin
      dat2 = dat2 * float(esi[indx[0]].exp / esi[indx[1]].exp)
      var2 = var2 * float(esi[indx[0]].exp / esi[indx[1]].exp)
  endif

  ;; Ratio
  rtio = dat1/dat2
  ;; Stats on the ratio
  djs_iterstat, rtio, sigrej=3.0, median=med_rtio, sigma=sig_rtio, maxiter=2

  ; Find all bad pixels
  bdpix = where(abs(rtio-med_rtio) GT 10.*sig_rtio, nbad)

  ; Average (not add)
  fimg = (dat1+dat2)/2.
  if arg_present(VAR) then var = (var1+var2)/2.


  ; Take minimum of bad pixels
  if nbad GT 0 then begin
      fimg[bdpix] = dat1[bdpix] < dat2[bdpix]
      if arg_present(VAR) then var[bdpix] = ( var1[bdpix] < var2[bdpix] )
  endif

  return, fimg

end
