;+
; NAME:
;   long_arc_fw
;
; PURPOSE:
;   Compute the full width of an arcline in an arc image corresponding to 
;   a given fraction of the peak flux
;
; CALLING SEQUENCE:
;    fw_frac=arc_fw(profile, xpeak, maxsep, frac)
; INPUTS:
;   profile      - 1-d arc spectrum
;   xpeak        - location of the arc line peak in pixels
;   maxsep       - separation to look for full width on either side
;   frac         - Fraction of peak flux to use for full width
;
; OUTPUTS: 
;   fw_frac      - Full width at fraction frac of peak flux
;
; OPTIONAL OUTPUTS: 
;
; COMMENTS:
; 
; EXAMPLES:
;
; BUGS:
;    
; PROCEDURES CALLED:
; 
;
; REVISION HISTORY:
;   27-May-2005 Written by J. Hennawi (UCB)
;-
;------------------------------------------------------------------------------
FUNCTION LONG_ARC_FWHM, profile, xpeak_in, maxsep

nx = n_elements(profile)
xarr = dindgen(nx)
i1 = round(xpeak_in-maxsep) > 0
i2 = round(xpeak_in+maxsep) < (nx-1L)
xpeak = long_find_nminima(-profile[i1:i2], xarr[i1:i2], nfind = 1, minsep = 3 $
                          , ypeak = ypeak, npeak = npeak, errcode = errcode $
                          , width = maxsep, sigma = sigma $
                          , /doplot, xplotfit = xfit, yplotfit = yfit)
fwhm_frac = 2.35482D*sigma[0]

RETURN, fwhm_frac
END
