;+
; NAME:
;   long_stripstd
;
; PURPOSE:
;
;   Create a trace FITS file from an existing std-.. file
;
; CALLING SEQUENCE:
;  long_reduce, planfile, /clobber, /NOZAP, /NOFLEX, /NOHELIO 
;
; INPUTS:
;  planfile  -- File created by long_plan which guides the reduction
;               process
;
; OPTIONAL INPUTS:
; /NOFLEX  -- Do not apply flexure correction [necessary if your setup
;             has not been calibrated.  Contact JH or JXP for help if
;             this is the case.]
;  HAND_FWHM -- Set the FWHM of the object profile to this value (in
;               pixels)
; /NOHELIO -- Do not correct to heliocentric velocities
; /NOZAP   -- Do not flag CRs
; REID_FILE=  -- File for cross-correlation [Replaces the default]
;
; OUTPUTS:
;  (1) Various calibration files
;  (2) One multi-extension FITS file in Science per exposure containing
;  the extracted data and processed images
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
; BUGS:
;   
; REVISION HISTORY:
;   12-Jun-2008  Written by JXP
;-  
;-----------------------------------------------------------------------------
PRO long_stripstd, stdfile, trcfile

  if N_params() LT 2 THEN BEGIN
      print,'Syntax - ' + $
        'long_stripstd, stdfile, trcfile  [v1.0]' 
      return
  endif 

  ;; Read in 
  stdstruct = xmrdfits(stdfile, 5, /silent)
  stdmax = max(stdstruct.PEAKFLUX, stdind)
  stdtrace = stdstruct[stdind]
  

  ;; Write
  print, 'long_stripstd: Writing ', trcfile
  dummy = replicate(1., 2,2)
  mwrfits, dummy, trcfile, /create ; 0
  mwrfits, dummy, trcfile          ; 1
  mwrfits, dummy, trcfile          ; 2
  mwrfits, dummy, trcfile          ; 3
  mwrfits, dummy, trcfile          ; 4
  mwrfits, stdtrace, trcfile       ; 5

  spawn, 'gzip -f '+trcfile

  return
end
;------------------------------------------------------------------------------
