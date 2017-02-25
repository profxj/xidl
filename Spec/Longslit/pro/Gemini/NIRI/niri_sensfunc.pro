;+
; NAME:
;   gnirs_sensfunc
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
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   05-June-2006  Written by J. Hennawi UC Berkeley
;-----------------------------------------------------------------------------

PRO niri_sensfunc, tellfiles, type, sensfuncfile, V = V, magfunc = magfunc $
                   , loglam = loglam, flux = flux, ivar = ivar  $
                   , CHECK = CHECK, OPT = OPT, COMBCHECK = COMBCHECK $
                   , objid = objid, arr_objid = objid_arr, IN_NPOLY = IN_NPOLY

IF NOT KEYWORD_SET(OPT) THEN BOX = 1
niri_fluxcal, tellfiles, loglam = loglam, flux = flux, ivar = ivar $
              , exptime = exptime_arr, box = box, CHECK = COMBCHECK $
              , objid = objid, arr_objid = objid_arr, IN_NPOLY = IN_NPOLY
ngrid = n_elements(flux)
exptime = total(exptime_arr)/double(n_elements(exptime_arr))
; Create telluric standard model spectrum with correct normalization
gnirs_telluric_std, type, loglam = loglam_std, flux = flux_std, V = V

flux_std = 1.0d17*flux_std      ; fluxes are in units of 1.0e-17
flux = flux/exptime
ivar = ivar*exptime^2
magfunc = fltarr(ngrid) - 100.0 ; big negative number
; find the min and max of the calibration spectrum
loglam_min_std = min(loglam_std)
loglam_max_std = max(loglam_std)
; find the min and max of the observed standard
ind = WHERE(loglam GT loglam_min_std AND $
            loglam LT loglam_max_std, ngood)
loglam1 = loglam[ind]
flux1   = flux[ind]
ivar1   = ivar[ind]
fluxlog  = dblarr(ngood)
magfunc1 = dblarr(ngood)
;   interpolate standard star spectrum onto observed wavelengths
flux_std_int = interpol(flux_std, loglam_std, loglam1)
pos_error = 1./sqrt((ivar1 > 0) + (ivar1 LT 0))
pos_mask = (flux1 GT pos_error/100.0) AND (ivar1 GT 0) $
  AND flux_std_int GT 0.0
pos = where(pos_mask, npos, COMPLEMENT = nopos, NCOMPLEMENT = nbad)
fluxlog[pos] = 2.5*alog10(flux1[pos] > (pos_error[pos]/100.0))
magfunc1[pos] = 2.5*alog10(flux_std_int[pos] > 1.0e-2) - fluxlog[pos]
;   sensfunc = 10.0^(0.4*magfunc)
magfunc[ind[pos]] = magfunc1[pos]
;   Now interpolate the magfunc for the masked pixels
IF nbad NE 0 THEN BEGIN
    magfunc_int = interpol(magfunc1[pos], loglam1[pos], loglam1[nopos])
    magfunc[ind[nopos]] = magfunc_int
ENDIF
; cap the magfunc so that sensfunc < 1.0e10
magfunc = magfunc <  15.0
IF keyword_set(CHECK) then x_splot, 10.0d^loglam/(1.0d4), 10.0D^(0.4D*magfunc) $
  , psym1 = 10, /blo 


;nocalib_inds = WHERE(loglam LE loglam_min_std OR loglam GE loglam_max_std, nno)
;IF nno NE 0 THEN magfunc[nocalib_inds] = 0.0

IF KEYWORD_SET(sensfuncfile) THEN BEGIN
    mwrfits, magfunc, sensfuncfile, /create
    mwrfits, loglam, sensfuncfile
ENDIF

RETURN
END
