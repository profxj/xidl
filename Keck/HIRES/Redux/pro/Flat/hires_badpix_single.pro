;+ 
; NAME:
; hires_badpix_single
;     Version 1.1
;
; PURPOSE:
;    Generates a bad pixel mask for the original chip
;
; CALLING SEQUENCE:
;   
;  hires_badpix_single, 
;
; INPUTS:
;   hires     -  HIRES structure
;   setup    -  Setup identifier 
;   [chip]   -  Blue (1), Green (2), Red (3), or multiple (Default:
;              [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per setup and side
;  (e.g. 'Flats/Flat/_B_01_T.fits.gz')
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite the final fits file
;   /USEBIAS - Use the bias frame in bias subtraction
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_mktflat, hires, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_subbias
;  xcombine
;  hires_delov
;
; REVISION HISTORY:
;   17-Apr-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_badpix_single, msk

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hires_mktflat, hires, setup, [chip], /CLOBBER, ' + $
        '/USEBIAS [v1.1]'
      return
  endif 

  ;; All of these pixels are in the rotated frame!        

  msk = fltarr(2048L,2048L)  ;; Unbinned frame

  ;; Bad columns
  msk[*,85:86] = 1. 
  msk[*,1127] = 1. 
  msk[*,2006:2007] = 1. 

  ;; Bleeding
  msk[1666:1704,1126:1200] = 1.
  msk[1108:1240,2005:*] = 1.
  
  ;; Ink spot
  msk[913:1008,946:1073] = 1.

  print, 'hires_badpix_single: All done!'
  return
end
