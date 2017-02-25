;+ 
; NAME:
; sdss_ewna1
;  Version 1.1
;
; PURPOSE:
;  Measures the EW of the MgII lines (rest values) in the SDSS
;
; CALLING SEQUENCE:
;  sdss_ewna1, wave, flux, sig, strct, ZABS=
;
; INPUTS:
;  wave  -- Wavelength array
;  flux  -- Flux array (normalized)
;  sig   -- Sigma array
;
; RETURNS:
;  strct -- MgII structure with MgII EW filled up
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  ZABS=  -- Absorption redshift (required)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   27-Feb-2004 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_ewna1, wave, flux, sig, strct, ZABS=zabs, TI2=ti2, ONLYTI2=onlyti2

  if  N_params() LT 4  then begin 
    print,'Syntax - ' + $
             'sdss_ewna1, wave, flux, sig, strct, ZABS= [v1.1]'
    return
  endif 

  if not keyword_set( ZABS ) then stop
  if keyword_set(ONLYTI2) then istrt = 3L else istrt = 0L

  ;; Structure
  strct={sdssmgiistrct}

  ;; Wavelengths
  rwave= [5891.5833d, 5897.5581]
  nlin = n_elements(rwave)
  obswv = (1.+zabs)*rwave

  ;; dwv
  dwv = wave - shift(wave,1)

  ;; Set box size for MgII
  boxw = (obswv[1]-obswv[0])/2.

  flux = flux < 1.1
  ;; Loop
  for ii=istrt,nlin-1 do begin
      strct.wrest[ii] = rwave[ii]

      ;; EW
      if ii LE 2 then begin
          ;; Edges
          mn = min(abs(wave-(obswv[ii]-boxw)),ilhs)
          mn = min(abs(wave-(obswv[ii]+boxw)),irhs)
          ;; EW
          strct.ew[ii] = total( (1.-flux[ilhs:irhs])*dwv[ilhs:irhs]) / (1.+zabs)
          strct.sigew[ii] = sqrt(total( (sig[ilhs:irhs]*dwv[ilhs:irhs])^2)) $
            / (1.+zabs)
      endif else begin  ;; Ti2
          ;; Edges
          mn = min(abs(wave-obswv[ii]),icen)
          idx = icen - 2 + lindgen(5)
          ;; EW
          strct.ew[ii] = total( (1.-flux[idx])*dwv[idx]) / (1.+zabs)
          strct.sigew[ii] = sqrt(total( (sig[idx]*dwv[idx])^2)) / (1.+zabs)
      endelse
  endfor


  return


end 


