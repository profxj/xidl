;+
; NAME:
;   tspec_fluxcombine
;
; PURPOSE:
;   Splices together the various orders into a final, 1D spectrum
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OPTIONAL OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   display     (T. Robishaw)
;   js_rdplot   (JDS)
;   traceset2xy (idlutils)
;   xy2traceset (idlutils)
;   splog       (idlutils)
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   05-June-2006  Written by J. Hennawi UC Berkeley
;-----------------------------------------------------------------------------
;;
;; tspec_fluxcombine, 'FSpec/flux-SBS1602+576.fits', outfile='FSpec/SBS1602+576F.fits'
PRO TSPEC_FLUXCOMBINE, infile, sensfuncfile, loglam = loglam $
                       , ONEDFLUX = flux, ONEDIVAR = IVAR $
                       , ONEDMASK = MASK $
                       , outfile = outfile $
                       , SIGREJ = SIGREJ1, plot=plot

IF KEYWORD_SET(SIGREJ1) THEN sigrej = sigrej1 ELSE SIGREJ = 3.0D
; Read in sensitivity function (only used for masking regions)
;magfunc  = mrdfits(sensfuncfile[0], 0)

influx = mrdfits(infile, 0, scihdr)
inivar = mrdfits(infile, 1)
inmask = mrdfits(infile, 2)
loglam = mrdfits(infile, 3)

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
      inivar[*, k] GT 0.       AND $
      influx[*, k] GT -1000.0    AND $
      influx[*, k] LT 1000.0     AND $
      finite(inivar[*, k]) EQ 1  AND $
      finite(influx[*, k]) EQ 1  ;AND $
      ;magfunc[*, k] GT -10.0   AND; $
      ;magfunc[*, k] LT  10.0    
    bad_ivar = WHERE(finite(inivar[*, k]) NE 1, nbad)
    bad_flux = WHERE(finite(influx[*, k]) NE 1, nbf)
    IF nbad NE 0 THEN inivar[bad_ivar, k] = 0.0D
    IF nbf  NE 0 THEN influx[bad_flux, k] = 0.0D
ENDFOR

;; Force the median of the orders to overlap
; med_flux = dblarr(norders)

wave = 10.0d^loglam * 1e4 ;; Ang


observ = strtrim(sxpar(scihdr, 'OBSERVAT'),2) ;; APO vs Palomar
case observ of 
   'APO': wvlim = [ [10010., 10595], [12000., 12350.], [99999,99999]]
   else: wvlim = [ [10010., 10595], [11640., 12275.], [14400., 14790]]
endcase


;; Tie order 0 to 1
ind_01 = WHERE(finalmask[*, 0] and finalmask[*,1] AND wave GT wvlim[0,0] AND wave LT wvlim[1,0])
djs_iterstat, influx[ind_01, 0]/influx[ind_01,1], sigrej = sigrej, median = scale_01
influx[*, 1] = influx[*, 1]*scale_01
inivar[*, 1] = inivar[*, 1]/scale_01^2

;; Tie order 1 to 2
ind_12 = WHERE(finalmask[*, 1] and finalmask[*,2] AND wave GT wvlim[0,1] AND wave LT wvlim[1,1])
djs_iterstat, influx[ind_12, 1]/influx[ind_12,2], sigrej = sigrej, median = scale_12
influx[*, 2] = influx[*, 2]*scale_12
inivar[*, 2] = inivar[*, 2]/scale_12^2

;; Tie order 2 to 3
if wvlim[0,2] LT 90000 then begin
   ind_23 = WHERE(finalmask[*, 2] and finalmask[*,3] AND wave GT wvlim[0,2] AND wave LT wvlim[1,2])
   djs_iterstat, influx[ind_23, 2]/influx[ind_23,3], sigrej = sigrej, median = scale_23
   influx[*, 3] = influx[*, 3]*scale_23
   inivar[*, 3] = inivar[*, 3]/scale_23^2
endif


;; Weights
sn2_weights = inivar
weights = sn2_weights*finalmask
wght_sum = total(weights, 2)
;;
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
    mwrfits, 1e4*10.0d^loglam, outfile ;; Ang
    print, 'tspec_fluxcombine: Final file is ', outfile
 ENDIF

if keyword_set(PLOT) then $
   x_specplot, outfile, infl=2, /bloc
   
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


