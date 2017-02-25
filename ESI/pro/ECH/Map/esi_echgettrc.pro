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
;   16-Sep-2004 Significant modification (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function esi_echgettrc, xval, yval, fin_fit, NPIX=npix, FITFIL=fitfil, $
                        OLDWAY=oldway

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'trc = esi_echgettrc(xval, yval, [fin_fit], NPIX=, ' + $
        'FITFIL=, /OLDWAY)  [v1.0]'
      return, -1
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
      diff = xhtrc[0] - xval
      trc0 = x_calcfit(findgen(npix), FITSTR=fin_fit[0])
      return, trc0-diff
  endif
      
  ;; Off to the right
  if xval GT xhtrc[nhtrc-1] then begin
      diff = xval - xhtrc[nhtrc-1] 
      trcN = x_calcfit(findgen(npix), FITSTR=fin_fit[nhtrc-1])
      return, trcN+diff
  endif

  ;; In the Good
  a = where(xval GT xhtrc AND xval LE shift(xhtrc,-1), na)
  if na NE 1 then begin
      print, 'esi_echgettrc: Shouldnt be here!!'
      stop
  endif
  
  indx = a[0]

  if keyword_set( OLDWAY ) then begin
      sep = xhtrc[indx+1] - xhtrc[indx]
      frac = (xval - xhtrc[indx])/sep
      
      ;; Grab the traces
      trc_low = x_calcfit(findgen(npix), FITSTR=fin_fit[indx])
      trc_hi = x_calcfit(findgen(npix), FITSTR=fin_fit[indx+1])
      
      ;; Final
      all_sep = trc_hi - trc_low
      fin = trc_low+frac*all_sep
  endif else begin
      ;; Set of 9
      ordr = indx / 9
      if (indx+1)/9 GT ordr OR $
        (indx-1)/9 LT ordr then print, 'esi_echgettrc: Warning!! Edge of order'
      ;; Parse x values
      xordr = xhtrc[lindgen(9)+9*ordr]

      ;; Fit
      if ordr EQ 0 then nord = 3 else nord = 5
      newfit = fin_fit[ordr*9]
      newfit.nord = nord
      svcoeff = dblarr(nord+1)
      for qq=0L,nord do begin
          fitstr = x_setfitstrct(NORD=1, HSIG=2.5, LSIG=2.5, /FLGREJ,$
                                 NITER=2, MAXREJ=3)
          coeff = dblarr(9)
          for j=0L,8 do coeff[j] = (*fin_fit[ordr*9+j].ffit)[qq]
          ;; Fit
          fit = x_fitrej(xordr, coeff, FITSTR=fitstr)
          
          ;; Interpolate
          svcoeff[qq] = x_calcfit(xval,FITSTR=fitstr)
      endfor
      *newfit.ffit = svcoeff

      ;; Trace
      fin = x_calcfit(findgen(npix),FITSTR=newfit)
  endelse
      

  return, fin
end
