;+ 
; NAME:
; x_setfitstrct   
;     Version 1.0
;
; PURPOSE:
;    Subtracts a scattered light model from the data
;    Outputs to the same image (usually OV/ov_esi####.fits)
;
; CALLING SEQUENCE:
;   
;  x_setfitstrct, esi, /DFLAT
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  Image with scattered light removed
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_setfitstrct, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

function x_setfitstrct, NITER=niter, MINPT=minpt, MAXREJ=maxrej, $
                        FUNC=func, NORD=nord, HSIG=hsig, LSIG=lsig, $
                        FLGREJ=flgrej

; Options
  if not keyword_set( HSIG ) then hsig = 3.
  if not keyword_set( LSIG ) then lsig = 3.
  if not keyword_set( FUNC ) then func = 'POLY'
  if not keyword_set( NORD ) then nord = 3L
  if not keyword_set( MAXREJ ) then maxrej = 100L
  if not keyword_set( NITER ) then niter = 1L
  if not keyword_set( MINPT ) then minpt = 1L

;  This routine defines the Fit structure

  tmp = {fitstrct}
  tmp.niter = niter
  tmp.minpt = minpt 
  tmp.maxrej = maxrej
  tmp.func = func
  tmp.nord = nord
  tmp.hsig = hsig
  tmp.lsig = lsig
  if keyword_set( FLGREJ ) then tmp.flg_rej = flgrej

  return, tmp

end
  
         
