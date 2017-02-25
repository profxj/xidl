;+ 
; NAME:
; long_skybkpts
;     Version 1.2
;
; PURPOSE:
;    Choose the placement of breakpoints for the sky spectrum
;
; CALLING SEQUENCE:
;   
;  long_skybkpts, esi, indx, /IFLAT, /REDDOV
;
; INPUTS:
;   esi     -  ESI structure
;   indx    -  Index values
;
; RETURNS:
;
; OUTPUTS:
;  Fully processed image
;
; OPTIONAL KEYWORDS:
;  /SVOV    - Save OV files
;  /REDOOV  - Redo OV subtraction
;  /CLOBBER - Clobber existing image
;  /SUBSCAT - Subtract scattered light from data image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   long_skybkpts, esi, [20L], /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;               Written by JFH
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION LONG_SKYBKPTS, piximg, bsp_min, nx, ny, goodpix

yarr = findgen(ny)## replicate(1.0, nx)
srt = sort(piximg[goodpix])
sky_pix = piximg[goodpix[srt]]
ypix = yarr[goodpix[srt]]
ymin = min(round(ypix))
ymax = max(round(ypix))
;samplmin = dindgen(ymax-ymin+1L)
;samplmax = dindgen(ymax-ymin+1L)
samplmin = fltarr(ymax-ymin+1L)
samplmax = fltarr(ymax-ymin+1L)
iflag = lonarr(ymax-ymin+1L) + 1L
FOR kk = ymin, ymax DO BEGIN
    smpix = where(ypix EQ kk, np)
    IF np GT 0 THEN BEGIN
        samplmin[kk-ymin] = min(sky_pix[smpix])
        samplmax[kk-ymin] = max(sky_pix[smpix])
    ENDIF ELSE iflag[kk-ymin] = 0
ENDFOR
samplmin = samplmin[where(iflag)]
samplmax = samplmax[where(iflag)]
nbkpt = n_elements(samplmax)
dsamp = shift(samplmin, -1)-samplmax
dsamp[nbkpt-1L] = dsamp[nbkpt-2L]
dsamp = djs_median(dsamp, width = 15, boundary = 'reflect')
dsamp = smooth(dsamp, 5)
;; if more than 60% of the pixels have dsamp < bsp_min than just use
;; uniform breakpoint spacing. Note that dsamp < 0 implies continuous
;; coverage. 
ilt = WHERE(dsamp LE bsp_min, nle)
IF double(nle) GE double(0.8*nbkpt) THEN BEGIN
    fullbkpt = bspline_bkpts(sky_pix, nord = 4, bkspace = bsp_min, /silent) 
    splog, 'Sampling of wavelength is nearly continuous. Using uniform spacing'
    splog, 'bsp = ', string(bsp_min, format = '(F5.2)')
    RETURN, fullbkpt
ENDIF
skybkpt_orig = samplmax + dsamp/2.0
skybkpt_orig = skybkpt_orig[sort(skybkpt_orig)]
;; compute the distance between breakpoints
dsamp2 = skybkpt_orig-shift(skybkpt_orig, 1)
dsamp2[0] = dsamp2[1]
igd = WHERE(dsamp2 GE bsp_min, nclose)
IF nclose GT 0 THEN skybkpt_orig = skybkpt_orig[igd]
skybkpt = skybkpt_orig
dsamp2 = skybkpt_orig-shift(skybkpt_orig, 1)
dsamp2[0] = dsamp2[1]
nbkpt=n_elements(skybkpt_orig)


;; Now loop over the bkpts and insert bkpts where we can add more
FOR kk = 1L, nbkpt-1L DO BEGIN
    dbkpt = skybkpt_orig[kk]-skybkpt_orig[kk-1L]
    ;; can we fit another bkpt
    dsamp_eff = dsamp2[kk] > bsp_min 
    IF dbkpt GT 2.0*dsamp_eff THEN BEGIN
        nsmp = floor(dbkpt/(dsamp_eff)) 
        bkpt_new = skybkpt_orig[kk-1L] + $
          (lindgen(nsmp-1L) + 1L)*dbkpt/double(nsmp)
        ibkpt = WHERE(skybkpt EQ skybkpt_orig[kk-1L])
        IF ibkpt EQ 0 THEN $
          skybkpt = [skybkpt[0], bkpt_new, skybkpt[ibkpt+1:*]] $
        ELSE IF ibkpt EQ n_elements(skybkpt)-2L THEN $
          skybkpt = $
          [skybkpt[0:ibkpt], bkpt_new, skybkpt[ibkpt+1]] $
        ELSE $
          skybkpt = $
          [skybkpt[0:ibkpt], bkpt_new, skybkpt[ibkpt+1:*]]
    ENDIF
ENDFOR

RETURN, skybkpt
END
