;+ 
; NAME:
; snr
;  V1.1
;
; PURPOSE:
;    Calculate the EW limit for an absorption feature 
;  for an nsig detection in an npix feature given a SNR per pix.
;
; CALLING SEQUENCE:
;   wlim, npix, snr, nsig, deltl
;
; INPUTS:
;  npix -- Number of pixels making up the feature (generally 3-4)
;  snr -- SNR per pixel
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
pro wlim, npix, snr, nsig, deltl
  a = deltl * sqrt(npix) * nsig / snr
  print, a
return
end
