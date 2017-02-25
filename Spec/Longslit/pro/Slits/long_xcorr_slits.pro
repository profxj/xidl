; NAME:
;   long_xcorr_slits
;
; PURPOSE:
;   Cross correlate two LRIS images for flexure compensation
;
; CALLING SEQUENCE:
;   xshift=lris_xcorr_slits(image1, tset_slits, maxlag=, /silent,/shift)
;
; INPUTS:
;   image1         - Reference frame 
;   tset_slits     - Slit set which will be shifted to align with image1
;
; OPTIONAL INPUTS: 
;   maxlag         - Maximum lag; default to 15 pixels
;   silent         - If set suppress reporting shift
;   shift          - If set the slit set will be returned shifted
; OUTPUTS:
;   xshift         - Number of pixels to shift slit set to align with IMAGE1,
;                    reported in fractional pixels
; COMMENTS:
;  
; EXAMPLES:
;   IDL> xshift=lris_xcorr_slits(image1,tset_slits)
;
; BUGS:
;    If the lags are not big enough this will return junk. Need to add
;    error checking. 
;    
; PROCEDURES CALLED:
;   find_npeaks()
;
; REVISION HISTORY:
;   May-2005  Written by D. Schlegel, J. Hennawi and S. Burles. 
;-
;------------------------------------------------------------------------------
FUNCTION long_xcorr_slits, image1, tset_slits, maxlag = maxlag $
                           , silent = silent, shift = shift $
                           , lags = lags, corr = corr
                           

if (size(tset_slits[0].COEFF, /n_dimen) EQ 1) then nslit = 1 $
ELSE nslit = (size(tset_slits[0].COEFF, /dimens))[1]

IF TAG_EXIST(tset_slits, 'XCORR_COEFF') THEN BEGIN
    if (size(tset_slits[0].XCORR_COEFF, /n_dimen) EQ 1) then nslit_corr = 1 $
    else nslit_corr = (size(tset_slits[0].XCORR_COEFF, /dimens))[1]
    xcorr_proto = $
      { func    :    tset_slits[0].FUNC, $
        xmin    :    tset_slits[0].XMIN, $
        xmax    :    tset_slits[0].XMAX, $
        coeff   :    dblarr(3, nslit_corr),   $
        dims    :    tset_slits[0].DIMS  $
      }
    xcorr_slits = replicate(xcorr_proto, 2)
    xcorr_slits[0].COEFF = tset_slits[0].XCORR_COEFF
    xcorr_slits[1].COEFF = tset_slits[1].XCORR_COEFF
ENDIF ELSE xcorr_slits = tset_slits
nx = tset_slits[0].dims[0]
ny = tset_slits[0].dims[1]

maskim = long_slits2mask(xcorr_slits, nslit = nslit1)
image2 = (maskim NE 0)

if (size(image1, /n_dimen) NE size(image2, /n_dimen) $
    OR total(size(image1, /dimens) NE size(image2, /dimens)) NE 0) then $
  message, 'ERROR: Dimensions of IMAGE1 and IMAGE2 must agree!'

if (NOT keyword_set(maxlag)) then maxlag = 15L
;lags = lindgen(2L * maxlag + 1L) - 15L
lags = lindgen(2L * maxlag + 1L) - maxlag   ;; KHRR changed this! -- Oct 6, 2014

igood = where(image2 GT 0)
djs_iterstat, image1[igood], mean = mean, median = median, sigma = sigma $
              , sigrej = 3.0
;;sky, image1*(image2 GT 0), skymode, skysig, /silent
;;thresh = skymode + 5.0d*skysig
thresh = median + 5.0d*sigma
corr = c2_correlate(image2, (image1 < thresh), lags)
;template = where(image2 GT 0)
;nx = (size(image2))[1]
;npix = n_elements(image2)
;template = template[where(template mod nx GT (-min(lags)) AND template mod nx LT (npix - max(lags)))]
;corr2 = lags*0.
;for i = 0, n_elements(lags)-1 do begin 
;    t = template + lags[i] 
;    corr2[i] = total(image1[t]) 
;endfor
xshift = long_find_nminima(-corr, lags, nfind = 1, width = 4, minsep = 2)
;xshift = find_npeaks(corr2, lags, nfind = 1)

if NOT (keyword_set(silent)) then $
  splog, 'Correlation gives shift of ', xshift

IF xshift EQ maxlag THEN begin
    print, 'WARNING: Cross correlation reached boundary. Increase lags'
    print, 'Or continue at your own risk'
    stop
    return, 0.
endif

; If shift is set shift both sets of coefficients
tset_slits_shift = long_shiftslits(tset_slits, xshift)

;slitmask_in  = long_slits2mask(tset_slits, nslit = nslit)
;slitmask_shf = long_slits2mask(tset_slits_shift)

;shiftmask = lonarr(nx, ny)
;FOR slitid = 1L, nslit DO BEGIN
    ;; these are the pixels that overlap between old slits and new slits
;    goodpix = WHERE(slitmask_in EQ slitid AND slitmask_shf EQ slitid)
;    shiftmask[goodpix] = 1L
;ENDFOR

IF KEYWORD_SET(SHIFT) THEN tset_slits = tset_slits_shift

RETURN, xshift
end
;------------------------------------------------------------------------------
