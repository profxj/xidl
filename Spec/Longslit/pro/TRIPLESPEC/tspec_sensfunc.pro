;+
; NAME:
;   tspec_sensfunc
;
; PURPOSE:
;   Use a standard star spectrum to determine the spectroscopic
;   response function and flux calibrate other spectra
;
; CALLING SEQUENCE:
;
; INPUTS:
;   scifile         - file containing object structure which has spectrum 
;                     for standard star
;
;   standard_name   - name of the standard star
;   
;   sensfuncfile    - File to write the sensitivity function out to
;
;
; OPTIONAL INPUTS:
;   OBJID           - object id of standar star in the object structure. Default
;                     is to the first object. 
;   irafdir         - the directory in which IRAF is installed
;                     (default set for Linux machines at Berkeley)
;   extinctfile     - file containing atmospheric extinction data as a
;                     function. Currently Mauna Kea and KPNO are supported
;                     and the extinction files will be read in automatically
;
; OPTIONAL OUTPUTS:
;   sensfunc        - structure containing sensitivity function
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;                   See README file in /apps2/iraf211/iraf/noao/lib/onedstds/
;                   for list of standard stars and the names of the
;                   associated files
;
; EXAMPLES:
;
; BUGS:
;                   Does not take into account atmospheric extinction!!!
;                   Leaves out first and last wavelength bins of
;                   sensitivity function
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

PRO tspec_sensfunc, tellfiles, type, sensfuncfile, V = V, magfunc = magfunc $
                    , loglam = loglam, flux = flux, ivar = ivar $
                    , itell = itell, OPT = OPT, CHECK = CHECK

IF NOT KEYWORD_SET(OPT) THEN BOX = 1
;IF NOT KEYWORD_SET(ITELL) THEN ITELL = 0L
scihdr = headfits(tellfiles[0])
; Rebin spectra onto common wavelength grid
tspec_fluxcal, tellfiles, loglam = loglam, flux = flux, ivar = ivar $
               , mask = mask, box = box;, ORDER_FIL=order_fil ;, /check
; Trying to combine the tellurics can result in systematic errors because 
; sometimes the negative trace is not on the slit. There should be plenty of 
; SNR in one exposure 
;flux = flux_arr[*, *, itell]
;ivar = ivar_arr[*, *, itell]

dims = size(flux, /dim)
ngrid = dims[0]
norders = dims[1]
;stop ;; The following line needs to be fixed -- JXP 24 April 2013
exptime = sxpar(scihdr, 'EXPTIME')
;exptime = 3.813
gnirs_telluric_std, type, loglam = loglam_std, flux = flux_std, V = V
flux_std = 1.0d17*flux_std      ; fluxes are in units of 1.0e-17
flux = flux/exptime
ivar = ivar*exptime^2

;; convert to microns
loglam_std = loglam_std - 4.

; these are the min max wavelengths of the orders (in microns)
;loglam_ord_min = [4.270, 4.146, 4.048, 3.970, 3.904, 3.840]
;loglam_ord_max = [4.410, 4.290, 4.190, 4.110, 4.040, 4.000]

observ = strtrim(sxpar(scihdr, 'OBSERVAT'),2) ;; APO vs Palomar
case observ of 
   'APO': begin
      loglam_ord_min = [-0.0905, -0.02477, 0.0541, 0.1506, 0.2748]
      loglam_ord_max = [0.02724, 0.09359, 0.1723, 0.2687, 0.3931]
   end
   else: begin ;; Palomar
      loglam_ord_min = [-0.0239, -0.01885, 0.0574, 0.150, 0.2731]
      loglam_ord_max = [0.02638, 0.09354, 0.1718, 0.2674, 0.3917]
   end
endcase

magfunc = fltarr(ngrid, norders) - 100.0 ; big negative number
; find the min and max of the calibration spectrum
loglam_min_std = min(loglam_std)
loglam_max_std = max(loglam_std)
; find the min and max of the observed standard
FOR k = 0L, norders-1L DO BEGIN
    ind_ord = WHERE(loglam GT loglam_ord_min[k] AND $
                    loglam LT loglam_ord_max[k], nord)
    loglam1 = loglam[ind_ord]
    flux1   = flux[ind_ord, k]
    ivar1   = ivar[ind_ord, k]
    fluxlog  = dblarr(nord)
    magfunc1 = dblarr(nord)
;   interpolate standard star spectrum onto observed wavelengths
    flux_std_int = interpol(flux_std, loglam_std, loglam1)
    pos_error = 1./sqrt((ivar1 > 0) + (ivar1 LT 0))
    pos_mask = (flux1 GT pos_error/10.0) AND (ivar1 GT 0) $
      AND flux_std_int GT 0.0
    pos = where(pos_mask, npos, COMPLEMENT = nopos, NCOMPLEMENT = nbad)
    fluxlog[pos] = 2.5*alog10(flux1[pos] > (pos_error[pos]/10))
    magfunc1[pos] = 2.5*alog10(flux_std_int[pos] > 1.0e-2) - fluxlog[pos]
;   sensfunc = 10.0^(0.4*magfunc)
    magfunc[ind_ord[pos], k] = magfunc1[pos] <  15.0D
;   Now interpolate the magfunc for the masked pixels
    IF nbad NE 0 THEN BEGIN
        magfunc_int = interpol(magfunc1[pos], loglam1[pos], loglam1[nopos])
        magfunc[ind_ord[nopos], k] = magfunc_int <  15.0D
    ENDIF
    IF KEYWORD_SET(CHECK) THEN  BEGIN 
        indgood = WHERE(ivar[*, k] GT 0.0, ngood)
        IF ngood GT 0 THEN BEGIN
            minx = min(loglam[indgood])
            maxx = max(loglam[indgood])
        ENDIF ELSE BEGIN
            minx = min(newloglam)
            maxx = max(newloglam)
        ENDELSE
        ind = where(magfunc[*, k] GT 0.0 AND magfunc[*, k] LT 15.0)
        djs_iterstat, 10.0d^(0.4D*magfunc[ind, k]) $
                      , median = median, sigma = sigma
        ymax = abs(median) + 6.0*sigma
        if keyword_set(CHECK) then $
           x_splot, 10.0d^loglam, 10.0D^(0.4D*magfunc[*, k]) $
                    , psym1 = 10, /blo $
                    , title = 'Sensitivity function for order # ' $
                    + strcompress(string(k), /rem) $
                    , xmnx = 10.0d^[minx, maxx], ymnx = [0.0, ymax]
    ENDIF
ENDFOR


;nocalib_inds = WHERE(loglam LE loglam_min_std OR loglam GE loglam_max_std, nno)
;IF nno NE 0 THEN magfunc[nocalib_inds] = 0.0

IF KEYWORD_SET(sensfuncfile) THEN BEGIN
   print, 'tspec_sensfunc: Writing ', sensfuncfile
   mwrfits, magfunc, sensfuncfile, /create
   mwrfits, loglam, sensfuncfile
ENDIF

RETURN
END
