;+ 
; NAME:
; esi_echgettrc   
;     Version 1.0
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;  esi_echgettrc, esi, /DFLAT
;
; INPUTS:
;   esi     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;  Image with gettrcered light removed
;
; OPTIONAL KEYWORDS:
;   DFLAT      - Use Dome flats where possible
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echgettrc, esi
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function esi_echgettrc, xval, yval, all_fit, NPIX=npix

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'trc = esi_echgettrc(xval, yval, all_fit)  [v1.0]'
      return, -1
  endif 
  
;  Optional Keywords
  if not keyword_set( NPIX ) then npix = 4096L
;  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'

;  Find traces

  nhtrc = n_elements(all_fit)

  xhtrc = fltarr(nhtrc)
  for q=0,nhtrc-1 do xhtrc[q] = x_calcfit(yval, FITSTR=all_fit[q])

  ;; Off to the left
  if xval LT xhtrc[0] then begin
      diff = xhtrc[0] - xval
      trc0 = x_calcfit(findgen(npix), FITSTR=all_fit[0])
      return, trc0-diff
  endif
      
  ;; Off to the right
  if xval GT xhtrc[nhtrc-1] then begin
      diff = xval - xhtrc[nhtrc-1] 
      trcN = x_calcfit(findgen(npix), FITSTR=all_fit[nhtrc-1])
      return, trcN+diff
  endif


  ;; In the Good
  a = where(xval GT xhtrc AND xval LE shift(xhtrc,-1), na)
  if na NE 1 then begin
      print, 'esi_echgettrc: Shouldnt be here!!'
      stop
  endif
  
  indx = a[0]

  sep = xhtrc[indx+1] - xhtrc[indx]
  frac = (xval - xhtrc[indx])/sep

  ;; Grab the traces
  trc_low = x_calcfit(findgen(npix), FITSTR=all_fit[indx])
  trc_hi = x_calcfit(findgen(npix), FITSTR=all_fit[indx+1])

  ;; Final
  all_sep = trc_hi - trc_low
  fin = trc_low+frac*all_sep

  return, fin
end
