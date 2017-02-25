;; Plots the output from tspec_fluxcal
PRO TSPEC_PLOT2D, flux_fil, ordr, ALL=all


;IF NOT KEYWORD_SET(ORDER_FIL) THEN stop
;; Generate order mask
;tset = xmrdfits(order_fil,1)
;omask = tspec_ordermask(tset, order_vec=ovec)

  IF KEYWORD_SET(ALL) THEN ordr = lindgen(5)
  IF NOT KEYWORD_SET(ORDR) THEN ordr = 0L
  nordr = n_elements(ordr)

  for qq=0L,nordr-1 do begin
     data = xmrdfits(flux_fil, 0)
     ivar = xmrdfits(flux_fil, 1)
     wave = 1e4 * (10.d^xmrdfits(flux_fil, 3))
     
     iordr = ordr[qq]
     fx = data[*,iordr]
     npix = n_elements(fx)
     sig = fltarr(npix)
     gd = where(ivar[*,iordr] GT 0)
     sig[gd] = 1./sqrt(ivar[gd,iordr])
     
     ;; Plot
     x_specplot, fx, sig, wave=wave, infl=4, /bloc
  endfor

return
end
