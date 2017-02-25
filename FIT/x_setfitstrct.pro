;+ 
; NAME:
; x_setfitstrct   
;     Version 1.1
;
; PURPOSE:
;    Initizlies a 1D FIT structure.
;
; CALLING SEQUENCE:
;   
;  fitstr = x_setfitstrct(NITER=, MINPT=, MAXREJ=, FUNC=, NORD=,
;  HSIG=, lSIG=, FLGREJ=)
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NITER -- Number of iterations for FIT  [default: 1L]
;  MINPT -- Number of points to keep in FIT [default: 1L]
;  MAXREJ -- Max Number of points to reject [default: 100L]
;  FUNC -- Max Number of points to reject [default: 'POLY']
;  NORD -- Order of the fit [default: 3L]
;  HSIG -- Upper sigma for rejection [default: 3.]
;  LSIG -- Lower sigma for rejection [default: 3.]
;  /FLGREJ -- Turn on rejection [default: NONE]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fitstr= x_setfitstrct(/FLGREJ)
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
  
         
