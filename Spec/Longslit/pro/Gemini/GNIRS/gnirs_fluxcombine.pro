PRO GNIRS_FLUXCOMBINE, infile, sensfuncfile, loglam = loglam $
                       , ONEDFLUX = flux, ONEDIVAR = IVAR $
                       , ONEDMASK = MASK $
                       , outfile = outfile $
                       , SIGREJ = SIGREJ1,HIZQSO=HIZQSO,NO_THRESH = NO_THRESH $
                       , NO_SENS_THRESH=NO_SENS_THRESH,NO_TIE=NO_TIE

IF KEYWORD_SET(SIGREJ1) THEN sigrej = sigrej1 ELSE SIGREJ = 3.0D
IF KEYWORD_SET(NO_THRESH) THEN MASK_THRESH = 1e30 ELSE MASK_THRESH=1000.0
IF KEYWORD_SET(NO_SENS_THRESH) THEN SENS_THRESH = 1e30 ELSE SENS_THRESH=10.0
; Read in sensitivity function (only used for masking regions)
magfunc  = xmrdfits(sensfuncfile[0], 0)

influx = xmrdfits(infile, 0, scihdr)
inivar = xmrdfits(infile, 1)
inmask = xmrdfits(infile, 2)
loglam = xmrdfits(infile, 3)
dims = size(influx, /dim)
ngrid = dims[0]
norders = dims[1]

finalmask = lonarr(ngrid, norders) + 1L
flux = dblarr(ngrid)
ivar = dblarr(ngrid)
mask = dblarr(ngrid)
;; Compute median SNR per order
;snr = dblarr(norders)
FOR k = 0, norders-1L DO BEGIN
    ;masklam = gnirs_ordermask(k)
    ;ind_ord = WHERE(inivar[*, k] GT 0.0 AND $
    ;                finite(inivar[*, k]) EQ 1 AND masklam, nind)
    ;djs_iterstat, influx[ind_ord, k]*sqrt(inivar[ind_ord, k]) $
    ;              , sigrej = sigrej, median = median_sn
    ;;snr[k] = median_sn
    finalmask[*, k] = inmask[*, k] AND $
      inivar[*, k] GT -1.0       AND $
      influx[*, k] GT -MASK_THRESH    AND $
      influx[*, k] LT  MASK_THRESH    AND $
      finite(inivar[*, k]) EQ 1  AND $
      finite(influx[*, k]) EQ 1  AND $
      magfunc[*, k] GT -SENS_THRESH   AND $
      magfunc[*, k] LT  SENS_THRESH    
    bad_ivar = WHERE(finite(inivar[*, k]) NE 1, nbad)
    bad_flux = WHERE(finite(influx[*, k]) NE 1, nbf)
    IF nbad NE 0 THEN inivar[bad_ivar, k] = 0.0D
    IF nbf  NE 0 THEN influx[bad_flux, k] = 0.0D
ENDFOR

;; Force the median of the orders to overlap
; med_flux = dblarr(norders)

wave = 10.0d^loglam
IF NOT KEYWORD_SET(NO_TIE) THEN BEGIN 
   ;; Tie order 4 to 3
   ind_3 = WHERE(finalmask[*, 0] AND wave GT 19400 AND wave LT 21200)
   ind_4 = WHERE(finalmask[*, 1] AND wave GT 17000 AND wave LT 18000)
   djs_iterstat, influx[ind_3, 0], sigrej = sigrej, median = median_3
   djs_iterstat, influx[ind_4, 1], sigrej = sigrej, median = median_4 
   scale_34 = (median_3/median_4)
   influx[*, 1] = influx[*, 1]*scale_34
   inivar[*, 1] = inivar[*, 1]/scale_34^2
   ;; Tie order 5 to 4
   ind_4 = WHERE(finalmask[*, 1] AND wave GT 14500 AND wave LT 16000)
   ind_5 = WHERE(finalmask[*, 2] AND wave GT 14300 AND wave LT 15000)
   djs_iterstat, influx[ind_4, 1], sigrej = sigrej, median = median_4
   djs_iterstat, influx[ind_5, 2], sigrej = sigrej, median = median_5 
   scale_45 = (median_4/median_5)
   influx[*, 2] = influx[*, 2]*scale_45
   inivar[*, 2] = inivar[*, 2]/scale_45^2
;; Tie order 6 to 5
   ind_5 = WHERE(finalmask[*, 2] AND wave GT 11900 AND wave LT 13000)
   ind_6 = WHERE(finalmask[*, 3] AND wave GT 11200 AND wave LT 12500)
   djs_iterstat, influx[ind_5, 2], sigrej = sigrej, median = median_5
   djs_iterstat, influx[ind_6, 3], sigrej = sigrej, median = median_6 
   scale_56 = (median_5/median_6)
   influx[*, 3] = influx[*, 3]*scale_56
   inivar[*, 3] = inivar[*, 3]/scale_56^2
;; Tie order 7 to 6
   IF KEYWORD_SET(HIZQSO) THEN BEGIN 
      ind_6 = WHERE(finalmask[*, 3] AND wave GT 10400 AND wave LT 11000)
      ind_7 = WHERE(finalmask[*, 4] AND wave GT  10400 AND wave LT 11000)
   ENDIF ELSE BEGIN
      ind_6 = WHERE(finalmask[*, 3] AND wave GT 10000 AND wave LT 11000)
      ind_7 = WHERE(finalmask[*, 4] AND wave GT  9750 AND wave LT 11000)
   ENDELSE
;; JFH changed for high-z QSO work. JFH 04/2017
   djs_iterstat, influx[ind_6, 3], sigrej = sigrej, median = median_6
   djs_iterstat, influx[ind_7, 4], sigrej = sigrej, median = median_7

   scale_67 = (median_6/median_7)
   influx[*, 4] = influx[*, 4]*scale_67
   inivar[*, 4] = inivar[*, 4]/scale_67^2
;; Tie order 8 to 7
   ind_7 = WHERE(finalmask[*, 4] AND wave GT  9000 AND wave LT  9600)
   ind_8 = WHERE(finalmask[*, 5] AND wave GT  8100 AND wave LT  9400)
   djs_iterstat, influx[ind_7, 4], sigrej = sigrej, median = median_7
   djs_iterstat, influx[ind_8, 5], sigrej = sigrej, median = median_8 
   scale_78 = (median_7/median_8)
   influx[*, 5] = influx[*, 5]*scale_78
   inivar[*, 5] = inivar[*, 5]/scale_78^2
ENDIF
;; Read in an archived bspline fit to the (S/N)^2 of each order 
;; obtained by fitting a fluxed telluric standard spectrum.
archive_fluxfile = getenv('LONGSLIT_DIR') + $
                   '/calib/flux/GNIRS/gnirs_flux_orders.fits'
;archive_fluxfile = '/Users/joe/REDUX/gnirs_redux/1249+0607/gnirs_flux_orders.fits'
sn2_archive = xmrdfits(archive_fluxfile, 0)
loglam_archive = xmrdfits(archive_fluxfile, 1)

sn2_weights = dblarr(ngrid, norders)
FOR iorder = 0L, norders-1L DO $
  sn2_weights[*, iorder] = $
  interpol(sn2_archive[*, iorder], loglam_archive, loglam)

weights = sn2_weights*finalmask
wght_sum = total(weights, 2)
flux = total(weights*influx, 2)/(wght_sum + (wght_sum EQ 0.0))
var =  double(inivar NE 0.0)/(inivar + (inivar EQ 0.0))
newvar =  total(weights^2*var, 2)/(wght_sum + (wght_sum EQ 0.0))^2
;; finalmask 1=good,0=bad
mask = total(finalmask, 2) GT 0.0
ivar = mask/(newvar + (newvar EQ 0.0))

IF keyword_set(OUTFILE) THEN BEGIN
    sxaddpar, scihdr, 'BITPIX', -32
    sxaddpar, scihdr, 'NAXIS', 1
    sxaddpar, scihdr, 'NAXIS1', n_elements(flux)
    sxdelpar, scihdr, 'NAXIS2'
    sxdelpar, scihdr, 'BZERO'
    sxdelpar, scihdr, 'BSCALE'
    mwrfits, flux, outfile, scihdr, /create
    giv = where(mask GT 0., ngiv)
    sig = 0*ivar - 1.0D
    sig[giv] = 1./sqrt(ivar[giv])
    mwrfits, sig, outfile
    mwrfits, 10.0d^loglam, outfile
    print, 'long_coadd: Final file is ', outfile
ENDIF
; IF KEYWORD_SET(OUTFILE) THEN BEGIN
;     mwrfits, flux, outfile, scihdr[10:*], /create
;     mwrfits, ivar, outfile
;     mwrfits, mask, outfile
;     mwrfits, loglam, outfile
;     mwrfits, newflux, outfile
;     mwrfits, newivar, outfile
;     mwrfits, newmask, outfile
; ENDIF

RETURN
END


