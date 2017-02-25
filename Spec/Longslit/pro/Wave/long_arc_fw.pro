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
FUNCTION LONG_ARC_FW, profile, xpeak, maxsep, frac

stop

maxsep = round(maxsep)
nx = n_elements(profile)
; find the peak value of the profile
i0 = round(xpeak)
ind = [i0-1, i0, i0+1]
ypeak = interpol(profile[ind], ind, xpeak)
yval = frac*ypeak
IF (i0 LE nx-1L) THEN BEGIN
    j2 = (where(profile[i0:i0+maxsep < (nx-1L)] LT yval))[0]
    IF j2 EQ -1 THEN iright = 0 $
    ELSE iright = interpol([i0 +j2-1, i0+j2 < (nx-1L)] $
                           , profile[i0+j2-1:i0+j2 < (nx-1L)], yval)
ENDIF ELSE iright = 0
IF (i0 GT 0) THEN BEGIN
    i1 = (reverse(where(profile[i0-maxsep > 0:i0] LT yval)))[0]
    IF i1 EQ -1 THEN ileft = 0 $
    ELSE BEGIN
        j1 = n_elements(profile[i0-maxsep > 0:i0])-1 -i1
        ileft = interpol([i0-j1 > 0, i0-j1+1], profile[i0-j1 > 0:i0-j1+1], yval)
    ENDELSE
ENDIF ELSE ileft = 0

IF ileft EQ 0 OR iright EQ 0 THEN fw_frac = 0.0 $
ELSE fw_frac = iright-ileft

RETURN, fw_frac
END
