;+ 
; NAME:
; snr
;  V1.1
;
; PURPOSE:
;    Calculate the signal-to-noise ratio needed to achieve an
;  nsig detection in an npix feature to a wlim EW limit.
;
; CALLING SEQUENCE:
;   snr, npix, wlim, nsig, deltl
;
; INPUTS:
;  npix -- Number of pixels making up the feature (generally 3-4)
;  wlim -- EW limit (Ang)
;  nsig -- Number of sigma limit
;  deltl -- Width of pix in Ang
;
; RETURNS:
;
; OUTPUTS:
;  Prints the SNR per pix needed to realize these limits
;
; OPTIONAL KEYWORDS:
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
;   Oct-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro snr, npix, wlim, nsig, deltl
  a = deltl * sqrt(npix) * nsig / wlim
  print, a
return
end
