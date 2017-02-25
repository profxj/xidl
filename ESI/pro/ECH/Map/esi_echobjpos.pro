;+ 
; NAME:
; esi_echobjpos   
;     Version 1.0
;
; PURPOSE:
;   Guess at the position of an object in other slits
;   given its position elsewhere.  Pass back traces as desired
;
; CALLING SEQUENCE:
;  esi_echobjpos, xval, yval, xall, TRC=
;
; INPUTS:
;   xval -- x-position in an order
;   yval -- y-position in an order
;
; RETURNS:
;
; OUTPUTS:
;  xall -- Position at yval in all orders
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  TRC=  -  Traces of the object in each order
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echobjpos, 1000., 2048., xall
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   16-Sep-2004 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echobjpos, xval, yval, xall, NPIX=npix, FITFIL=fitfil, TRC=trc

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echobjpos, xval, yval, xall, NPIX=, FITFIL=, TRC=  [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( NPIX ) then npix = 4096L
  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'

;  Find traces
  if not keyword_set( FIN_FIT ) then restore, fitfil

  nhtrc = n_elements(fin_fit)

  xhtrc = fltarr(nhtrc)
  for q=0,nhtrc-1 do xhtrc[q] = x_calcfit(yval, FITSTR=fin_fit[q])

  ;;; EXTREME EDGES ;;;
  ;; Off to the left
  if xval LT xhtrc[0] then begin
      stop
      diff = xhtrc[0] - xval
      trc0 = x_calcfit(findgen(npix), FITSTR=fin_fit[0])
      return
  endif
      
  ;; Off to the right
  if xval GT xhtrc[nhtrc-1] then begin
      stop
      diff = xval - xhtrc[nhtrc-1] 
      trcN = x_calcfit(findgen(npix), FITSTR=fin_fit[nhtrc-1])
      return
  endif

  ;; In the Good
  a = where(xval GT xhtrc AND xval LE shift(xhtrc,-1), na)
  if na NE 1 then begin
      print, 'esi_echgettrc: Shouldnt be here!!'
      stop
  endif
  
  indx = a[0]
  ordr = indx / 9
  if (indx+1)/9 GT ordr OR (indx-1)/9 LT ordr then stop

  pos = indx - ordr*9
  sep = xhtrc[indx+1] - xhtrc[indx]
  frac = (xval - xhtrc[indx])/sep
      
  xall = fltarr(10)
  if arg_present(TRC) then TRC = fltarr(npix,10L)
  for qq=0L,9 do begin
      ;; Sep and frac
      ii = qq*9 + pos
      sep = xhtrc[ii+1] - xhtrc[ii]
      xall[qq] = frac*sep + xhtrc[ii]

      ;; TRC
      if keyword_set(TRC) then trc[*,qq] = esi_echgettrc(xall[qq],yval,$
                                                         fin_fit)
  endfor
      
  return
end
